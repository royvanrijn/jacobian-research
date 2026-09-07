#!/usr/bin/env python3
"""Align finite q12 pencils with the exact quadratic descendant.

The exact connected pencil was compiled over ``QQ(a,b)`` with ``theta=a+b``.
Every displayed moving-equation coefficient is in ``QQ(theta^2)=QQ(ab)``.
Write ``omega=4*a*b``; then ``omega^2=16*q1*q2``.  At each collected prime
this worker fixes the unique local generator scale

    omega -> kappa * (2*r + ell)

and the induced symmetric-square base-PGL2 gauge that makes all 63 exact
moving-equation coefficients agree with the independently compiled local
pencil.
"""

import argparse
import hashlib
import itertools
import json
import re
import sys
from pathlib import Path


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
BASE_PRIMES = (19, 61, 67, 83, 89, 103, 131)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--exact-pencil", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json",
)
parser.add_argument(
    "--operands", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json",
)
for prime in BASE_PRIMES:
    infix = "common-" if prime == 19 else ""
    parser.add_argument(
        f"--p{prime}", type=Path,
        default=RESULTS / f"q80-third-q12-um2-p{prime}-{infix}resolved-pencil.json",
    )
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-third-q12-um2-exact-quadratic-pencils-aligned.json",
)
parser.add_argument(
    "--extra-prime",
    action="append",
    nargs=2,
    metavar=("PRIME", "RESOLVED_PENCIL"),
    default=[],
    help="append an independently compiled inert good prime",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
args.exact_pencil = args.exact_pencil.resolve()
args.operands = args.operands.resolve()
args.output = args.output.resolve()
extra_paths = {}
for prime_text, path_text in args.extra_prime:
    prime = int(prime_text)
    if prime in BASE_PRIMES or prime in extra_paths:
        raise ValueError(f"duplicate extra prime {prime}")
    extra_paths[prime] = Path(path_text).resolve()
PRIMES = BASE_PRIMES + tuple(sorted(extra_paths))
local_paths = {prime: getattr(args, f"p{prime}").resolve() for prime in BASE_PRIMES}
local_paths.update(extra_paths)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


exact = json.loads(args.exact_pencil.read_text())
operands = json.loads(args.operands.read_text())
locals_by_prime = {
    prime: json.loads(path.read_text()) for prime, path in local_paths.items()
}
if exact.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_RESOLVED_PENCIL":
    raise ValueError("exact resolved pencil is not certified")
if operands.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL":
    raise ValueError("exact closure operands are not certified")
for prime, payload in locals_by_prime.items():
    allowed_statuses = {"PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER"}
    if prime == 19:
        allowed_statuses.add("PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC")
    if payload.get("status") not in allowed_statuses:
        raise ValueError(f"p={prime}: local resolved pencil is not certified")
    if payload["specialization"]["prime"] != prime:
        raise ValueError(f"p={prime}: specialization mismatch")


def parse_rational_string(value):
    if "/" in value:
        numerator, denominator = value.split("/")
        return int(numerator), int(denominator)
    return int(value), 1


exact_atoms = {}
for t_degree, w_degree, x_degree, encoded in exact["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    if len(encoded) != 1:
        raise ArithmeticError("exact coefficient is not one theta expression")
    match = re.fullmatch(r"(.+)\*theta\^2 ([+-]) (.+)", encoded[0])
    if match is None:
        raise ArithmeticError("exact pencil does not descend coefficientwise to QQ(theta^2)")
    alpha = parse_rational_string(match.group(1))
    beta = parse_rational_string(match.group(3))
    if match.group(2) == "-":
        beta = (-beta[0], beta[1])
    exact_atoms[(t_degree, w_degree, x_degree)] = (alpha, beta)
if len(exact_atoms) != 63:
    raise ArithmeticError("exact moving-equation support changed")


def rational_record_mod(record, prime):
    return (
        int(record["numerator"]) % prime
        * pow(int(record["denominator"]) % prime, -1, prime)
        % prime
    )


def rational_pair_mod(pair, prime):
    return pair[0] % prime * pow(pair[1] % prime, -1, prime) % prime


def field_operations(prime, linear, constant):
    """Arithmetic in F_p[r]/(r^2+linear*r+constant), as (1,r)."""

    zero = (0, 0)
    one = (1, 0)

    def add(left, right):
        return ((left[0] + right[0]) % prime, (left[1] + right[1]) % prime)

    def neg(value):
        return ((-value[0]) % prime, (-value[1]) % prime)

    def multiply(left, right):
        return (
            (left[0] * right[0] - constant * left[1] * right[1]) % prime,
            (
                left[0] * right[1]
                + left[1] * right[0]
                - linear * left[1] * right[1]
            ) % prime,
        )

    def inverse(value):
        norm = (
            value[0] ** 2
            - linear * value[0] * value[1]
            + constant * value[1] ** 2
        ) % prime
        inverse_norm = pow(norm, -1, prime)
        return (
            (value[0] - linear * value[1]) * inverse_norm % prime,
            -value[1] * inverse_norm % prime,
        )

    def total(values):
        result = zero
        for value in values:
            result = add(result, value)
        return result

    return zero, one, add, neg, multiply, inverse, total


def matrix_multiply(left, right, multiply, total):
    return [
        [
            total(multiply(left[row][index], right[index][column]) for index in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def matrix_inverse(matrix, zero, one, add, neg, multiply, inverse):
    size = len(matrix)
    augmented = [
        list(matrix[row]) + [one if row == column else zero for column in range(size)]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column] != zero),
            None,
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_inverse = inverse(augmented[column][column])
        augmented[column] = [multiply(pivot_inverse, value) for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                add(value, neg(multiply(factor, pivot_value)))
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [row[size:] for row in augmented]


q1_record = operands["biquadratic_field"]["q1"]
q2_record = operands["biquadratic_field"]["q2"]
alignments = {}
for prime in PRIMES:
    local = locals_by_prime[prime]
    modulus = local["specialization"]["extension_modulus"]
    match = re.fullmatch(r"r\^2\s*\+\s*(\d+)\*r\s*\+\s*(\d+)", modulus)
    if match is None:
        raise ValueError(f"p={prime}: cannot parse local quadratic modulus")
    linear, constant = map(int, match.groups())
    local_discriminant = (linear**2 - 4 * constant) % prime
    if pow(local_discriminant, (prime - 1) // 2, prime) != prime - 1:
        raise ArithmeticError(f"p={prime}: local quadratic discriminant is not a nonsquare")
    q1 = rational_record_mod(q1_record, prime)
    q2 = rational_record_mod(q2_record, prime)
    global_discriminant = 16 * q1 * q2 % prime
    if pow(global_discriminant, (prime - 1) // 2, prime) != prime - 1:
        raise ArithmeticError(f"p={prime}: global QQ(ab) field does not remain quadratic")
    scale_square = global_discriminant * pow(local_discriminant, -1, prime) % prime
    scale_candidates = [
        value for value in range(prime) if value * value % prime == scale_square
    ]
    if len(scale_candidates) != 2:
        raise ArithmeticError(f"p={prime}: local generator scale does not have two signs")

    zero, one, add, neg, multiply, inverse, total = field_operations(
        prime, linear, constant
    )
    local_terms = {
        (w_degree, x_degree, t_degree): tuple(coefficient)
        for t_degree, w_degree, x_degree, coefficient
        in local["moving_equation"]["terms_T_W_x_coefficient_1_r"]
    }
    labels = sorted({(w_degree, x_degree) for _, w_degree, x_degree in exact_atoms})
    if len(local_terms) != 63 or len(labels) != 21:
        raise ArithmeticError(f"p={prime}: local moving-equation support changed")

    hits = []
    inverse_two = pow(2, -1, prime)
    for scale in scale_candidates:
        exact_local = {}
        for label, (alpha_record, beta_record) in exact_atoms.items():
            t_degree, w_degree, x_degree = label
            alpha = rational_pair_mod(alpha_record, prime)
            beta = rational_pair_mod(beta_record, prime)
            # theta^2=q1+q2+omega/2 and omega=scale*(2r+linear).
            exact_local[(w_degree, x_degree, t_degree)] = (
                (
                    beta
                    + alpha
                    * (q1 + q2 + scale * linear * inverse_two)
                ) % prime,
                alpha * scale % prime,
            )
        exact_vectors = [
            [exact_local[(w_degree, x_degree, degree)] for w_degree, x_degree in labels]
            for degree in range(3)
        ]
        local_vectors = [
            [local_terms[(w_degree, x_degree, degree)] for w_degree, x_degree in labels]
            for degree in range(3)
        ]

        transformation = None
        for indices in itertools.combinations(range(len(labels)), 3):
            exact_minor = [
                [exact_vectors[row][column] for column in indices]
                for row in range(3)
            ]
            exact_minor_inverse = matrix_inverse(
                exact_minor, zero, one, add, neg, multiply, inverse
            )
            if exact_minor_inverse is None:
                continue
            local_minor = [
                [local_vectors[row][column] for column in indices]
                for row in range(3)
            ]
            transformation = matrix_multiply(
                local_minor, exact_minor_inverse, multiply, total
            )
            break
        if transformation is None:
            raise ArithmeticError(f"p={prime}: exact coefficient vectors have rank below three")
        replay = matrix_multiply(transformation, exact_vectors, multiply, total)
        if replay != local_vectors:
            continue

        # Coefficients of a binary quadratic have discriminant x1^2-4*x0*x2.
        # Its similitude group in projective dimension two is the symmetric
        # square image of PGL2, so this certifies the base-gauge interpretation.
        minus_two = ((-2) % prime, 0)
        conic_form = [
            [zero, zero, minus_two],
            [zero, one, zero],
            [minus_two, zero, zero],
        ]
        transpose = [list(column) for column in zip(*transformation)]
        transformed_form = matrix_multiply(
            matrix_multiply(transpose, conic_form, multiply, total),
            transformation,
            multiply,
            total,
        )
        similitude = transformed_form[1][1]
        expected_form = [
            [multiply(similitude, value) for value in row] for row in conic_form
        ]
        if transformed_form != expected_form or similitude == zero:
            raise ArithmeticError(f"p={prime}: coefficient transformation is not induced by PGL2")
        # Recover one explicit 2x2 representative.  For
        # V_exact=(a+b*V_local)/(c+d*V_local), clearing the square denominator
        # acts on [x0,x1,x2] by the displayed symmetric-square matrix.  Fix
        # c=1; all seven controls have a nonzero upper-left entry.
        quadratic_scalar = transformation[0][0]
        if quadratic_scalar == zero:
            raise ArithmeticError(f"p={prime}: selected PGL2 recovery chart is singular")
        divide_by_scalar = inverse(quadratic_scalar)
        a_entry = multiply(transformation[0][1], divide_by_scalar)
        two_scalar = multiply((2 % prime, 0), quadratic_scalar)
        d_entry = multiply(transformation[1][0], inverse(two_scalar))
        b_entry = add(
            multiply(transformation[1][1], divide_by_scalar),
            neg(multiply(a_entry, d_entry)),
        )
        c_entry = one
        reconstructed = [
            [
                multiply(quadratic_scalar, multiply(c_entry, c_entry)),
                multiply(quadratic_scalar, multiply(a_entry, c_entry)),
                multiply(quadratic_scalar, multiply(a_entry, a_entry)),
            ],
            [
                multiply(quadratic_scalar, multiply((2 % prime, 0), multiply(c_entry, d_entry))),
                multiply(
                    quadratic_scalar,
                    add(multiply(a_entry, d_entry), multiply(b_entry, c_entry)),
                ),
                multiply(quadratic_scalar, multiply((2 % prime, 0), multiply(a_entry, b_entry))),
            ],
            [
                multiply(quadratic_scalar, multiply(d_entry, d_entry)),
                multiply(quadratic_scalar, multiply(b_entry, d_entry)),
                multiply(quadratic_scalar, multiply(b_entry, b_entry)),
            ],
        ]
        if reconstructed != transformation:
            raise ArithmeticError(f"p={prime}: explicit PGL2 representative does not replay")
        determinant = add(
            multiply(a_entry, d_entry), neg(multiply(b_entry, c_entry))
        )
        if determinant == zero:
            raise ArithmeticError(f"p={prime}: recovered base transformation is singular")
        hits.append(
            (
                scale,
                transformation,
                similitude,
                quadratic_scalar,
                [[a_entry, b_entry], [c_entry, d_entry]],
                determinant,
            )
        )

    if len(hits) != 1:
        raise ArithmeticError(
            f"p={prime}: expected one exact generator/base-gauge alignment, got {len(hits)}"
        )
    (
        scale,
        transformation,
        similitude,
        quadratic_scalar,
        pgl2_matrix,
        pgl2_determinant,
    ) = hits[0]
    alignments[str(prime)] = {
        "local_extension_modulus": modulus,
        "local_anti_invariant_generator": f"2*r+{linear}",
        "local_generator_discriminant": local_discriminant,
        "global_omega_square_mod_p": global_discriminant,
        "unique_omega_to_local_generator_scale": scale,
        "omega_reduction": f"omega -> {scale}*(2*r+{linear})",
        "quadratic_V_coefficient_transformation_exact_to_local": [
            [[int(value[0]), int(value[1])] for value in row]
            for row in transformation
        ],
        "discriminant_conic_similitude": [
            int(similitude[0]), int(similitude[1])
        ],
        "base_coordinate_formula": (
            "V_exact=(a+b*V_local)/(c+d*V_local); the degree-two equation is multiplied by "
            "quadratic_scalar*(c+d*V_local)^2"
        ),
        "base_PGL2_matrix_a_b_c_d": [
            [[int(value[0]), int(value[1])] for value in row]
            for row in pgl2_matrix
        ],
        "base_PGL2_determinant": [
            int(pgl2_determinant[0]), int(pgl2_determinant[1])
        ],
        "quadratic_equation_scalar": [
            int(quadratic_scalar[0]), int(quadratic_scalar[1])
        ],
        "all_63_moving_coefficients_replayed": True,
        "base_transformation_is_symmetric_square_PGL2": True,
    }


output = {
    "schema": "elkies-k3.q80-third-q12-exact-quadratic-pencil-prime-alignment.v1",
    "status": "PASS_EXACT_QQ_THIRD_Q12_QUADRATIC_PENCIL_DESCENT_AND_LOCAL_ALIGNMENT",
    "specialization": {"u": "-2", "primes": list(PRIMES)},
    "exact_quadratic_descent": {
        "biquadratic_field": "QQ(a,b), a^2=q1, b^2=q2",
        "compiler_primitive_element": "theta=a+b",
        "pencil_coefficient_field": "QQ(theta^2)=QQ(a*b)",
        "eta": "theta^2=q1+q2+2*a*b",
        "omega": "2*eta-2*(q1+q2)=4*a*b",
        "omega_square": "16*q1*q2",
        "minimal_polynomial_eta": "eta^2-2*(q1+q2)*eta+(q1-q2)^2",
        "exact_moving_coefficients_in_quadratic_subfield": 63,
    },
    "local_alignments": alignments,
    "inputs": {
        "exact_pencil": {
            "path": str(args.exact_pencil.relative_to(ROOT)),
            "sha256": sha256(args.exact_pencil),
        },
        "operands": {
            "path": str(args.operands.relative_to(ROOT)),
            "sha256": sha256(args.operands),
        },
        "local_pencils": {
            str(prime): {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
            }
            for prime, path in local_paths.items()
        },
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "the exact 63-term connected pencil descends from QQ(a,b) to QQ(a*b)",
            f"the same global quadratic field remains inert at all {len(PRIMES)} collected primes",
            "a unique signed local generator scale is fixed at every prime by literal coefficient replay",
            "a base-PGL2 gauge induced on quadratic V coefficients aligns every local pencil with the exact pencil",
        ],
        "not_proved": [
            "base-gauge transport of the local Jacobian and birational-map coefficients",
            "CRT or rational reconstruction of a characteristic-zero child",
            "a characteristic-zero minimal model, fibre marking, or Mordell--Weil rank",
        ],
    },
    "reproduce": " ".join(
        [
            "python3",
            "elkies-k3/scripts/align_q80_third_q12_exact_quadratic_pencil_primes.py",
            "--exact-pencil", str(args.exact_pencil.relative_to(ROOT)),
            "--operands", str(args.operands.relative_to(ROOT)),
        ]
        + [
            item
            for prime in BASE_PRIMES
            for item in (f"--p{prime}", str(local_paths[prime].relative_to(ROOT)))
        ]
        + [
            item
            for prime in sorted(extra_paths)
            for item in (
                "--extra-prime",
                str(prime),
                str(local_paths[prime].relative_to(ROOT)),
            )
        ]
        + ["--output", str(args.output.relative_to(ROOT))]
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"quadratic pencil alignment artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12QUADRATICALIGN|field=QQ(ab)|omega=4ab|primes={','.join(map(str, PRIMES))}|"
    f"coefficients=63x{len(PRIMES)}|local_generators=unique|base_gauges=PGL2|"
    "status=PASS_EXACT_QQ_THIRD_Q12_QUADRATIC_PENCIL_DESCENT_AND_LOCAL_ALIGNMENT"
)
