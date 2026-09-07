#!/usr/bin/env sage -python
"""Run the pinned p=19 factor lift with resumable full-target checkpoints.

The mathematical worker stays byte-for-byte pinned.  This launcher checks its
hash, injects only checkpoint serialization around the two Hensel loops, and
executes it with the canonical ``__file__`` so the final certificate remains a
literal output of the pinned worker.  Checkpoints store canonical rational
functions modulo the *final* target precision; they never coerce a capped
lower-precision p-adic fraction field into a larger one.
"""

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CANONICAL = ROOT / "elkies-k3/scripts/lift_q80_third_q12_discriminant_factors_p19_adic.sage"
CANONICAL_SHA256 = "9b81a3dc5f599dc7b0c179f8390536f59d7a09771ac6caa495bb7321a9c4f4c4"


def replace_once(source, old, new, label):
    if source.count(old) != 1:
        raise RuntimeError(f"checkpoint injection anchor {label!r} is not unique")
    return source.replace(old, new, 1)


source = CANONICAL.read_text()
actual_sha256 = hashlib.sha256(source.encode()).hexdigest()
if actual_sha256 != CANONICAL_SHA256:
    raise RuntimeError(
        "canonical p-adic factor worker changed: "
        f"expected {CANONICAL_SHA256}, got {actual_sha256}"
    )

source = replace_once(
    source,
    'parser.add_argument("--check", action="store_true")\n',
    '''parser.add_argument("--check", action="store_true")
parser.add_argument("--checkpoint", type=Path)
parser.add_argument("--resume-checkpoint", type=Path)
parser.add_argument("--checkpoint-every", type=int, default=16)
parser.add_argument(
    "--pointwise-newton",
    action="store_true",
    help="use full-precision Newton solves at good U-values and interpolate",
)
parser.add_argument(
    "--stop-after",
    help="testing/controlled-stop gate in the form factors:N or root:N",
)
''',
    "arguments",
)
source = replace_once(
    source,
    'args.output = args.output.resolve()\nif args.digits < 2:\n',
    '''args.output = args.output.resolve()
if args.checkpoint is not None:
    args.checkpoint = args.checkpoint.resolve()
if args.resume_checkpoint is not None:
    args.resume_checkpoint = args.resume_checkpoint.resolve()
    if args.checkpoint is None:
        args.checkpoint = args.resume_checkpoint
if args.checkpoint_every < 1:
    raise ValueError("checkpoint interval must be positive")
checkpoint_stop = None
if args.stop_after is not None:
    stop_parts = args.stop_after.split(":")
    if len(stop_parts) != 2 or stop_parts[0] not in {"factors", "root"}:
        raise ValueError("--stop-after must be factors:N or root:N")
    checkpoint_stop = (stop_parts[0], int(stop_parts[1]))
    if checkpoint_stop[1] < 1:
        raise ValueError("checkpoint stop valuation must be positive")
if args.digits < 2:
''',
    "argument validation",
)

checkpoint_helpers = r'''

checkpoint_payload = None
if args.resume_checkpoint is not None:
    checkpoint_payload = json.loads(args.resume_checkpoint.read_text())
    if checkpoint_payload.get("schema") != "elkies-k3.q80-third-q12-p19-hensel-checkpoint.v1":
        raise ValueError("unrecognized p-adic checkpoint schema")
    if checkpoint_payload["canonical_worker_sha256"] != sha256(Path(__file__).resolve()):
        raise ValueError("checkpoint was not produced by the pinned canonical worker")
    if checkpoint_payload["input_sha256"] != sha256(args.input):
        raise ValueError("checkpoint source pencil differs from the requested input")
    if checkpoint_payload["target_digits"] != args.digits:
        raise ValueError("checkpoint target precision differs from --digits")

checkpoint_modulus = prime**working_precision


def checkpoint_constant_coordinates(value):
    coordinates = padic_field(value).polynomial().list()
    coordinates += [padic.zero()] * (2 - len(coordinates))
    return [int(coordinates[index].lift()) % checkpoint_modulus for index in range(2)]


def checkpoint_rational_record(value):
    value = padic_u_field(value)
    numerator = padic_u_ring(value.numerator())
    denominator = padic_u_ring(value.denominator())
    leading = denominator.leading_coefficient()
    numerator /= leading
    denominator /= leading
    return {
        "numerator": [checkpoint_constant_coordinates(value) for value in numerator.list()],
        "denominator": [checkpoint_constant_coordinates(value) for value in denominator.list()],
    }


def checkpoint_w_record(value):
    return [checkpoint_rational_record(coefficient) for coefficient in value.list()]


def load_checkpoint_constant(coordinates):
    return padic_field(padic(int(coordinates[0]))) + padic_field(
        padic(int(coordinates[1]))
    ) * omega


def load_checkpoint_rational(record):
    numerator = padic_u_ring(
        [load_checkpoint_constant(coordinates) for coordinates in record["numerator"]]
    )
    denominator = padic_u_ring(
        [load_checkpoint_constant(coordinates) for coordinates in record["denominator"]]
    )
    return padic_u_field(numerator) / padic_u_field(denominator)


def load_checkpoint_w(record):
    return padic_w_ring([load_checkpoint_rational(value) for value in record])


if checkpoint_payload is not None:
    L = load_checkpoint_w(checkpoint_payload["L"])
    Q = load_checkpoint_w(checkpoint_payload["Q"])
    D = load_checkpoint_w(checkpoint_payload["D"])


def write_checkpoint(stage, valuation):
    if args.checkpoint is None:
        return
    forced_stop = checkpoint_stop == (stage, valuation)
    if valuation < args.digits and valuation % args.checkpoint_every and not forced_stop:
        return
    record = {
        "schema": "elkies-k3.q80-third-q12-p19-hensel-checkpoint.v1",
        "canonical_worker_sha256": sha256(Path(__file__).resolve()),
        "checkpoint_runner_sha256": hashlib.sha256(
            Path(__checkpoint_runner_path__).read_bytes()
        ).hexdigest(),
        "input_path": str(args.input.relative_to(ROOT)),
        "input_sha256": sha256(args.input),
        "target_digits": args.digits,
        "working_precision": working_precision,
        "stage": stage,
        "valuation": valuation,
        "factor_valuation_history_prefix": valuation_history[:-1],
        "root_valuation_history_prefix": (
            root_valuation_history[:-1] if stage == "root" else []
        ),
        "L": checkpoint_w_record(L),
        "Q": checkpoint_w_record(Q),
        "D": checkpoint_w_record(D),
        "repeated_root": (
            checkpoint_w_record(repeated_root) if stage == "root" else None
        ),
    }
    serialized_checkpoint = json.dumps(record, indent=2, sort_keys=True) + "\n"
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.checkpoint.with_suffix(args.checkpoint.suffix + ".tmp")
    temporary.write_text(serialized_checkpoint)
    temporary.replace(args.checkpoint)
    print(
        f"Q80THIRDQ12PADICCHECKPOINT|stage={stage}|valuation={valuation}|"
        f"path={args.checkpoint}",
        flush=True,
    )
    if forced_stop:
        raise SystemExit(0)


checkpoint_local_w_ring = PolynomialRing(padic_field, "W_checkpoint")
W_checkpoint = checkpoint_local_w_ring.gen()


def evaluate_rational_at(value, point):
    value = padic_u_field(value)
    denominator = value.denominator()(point)
    if not denominator or denominator.valuation() != 0:
        raise ZeroDivisionError
    return value.numerator()(point) / denominator


def evaluate_w_at(value, point):
    return checkpoint_local_w_ring(
        [evaluate_rational_at(coefficient, point) for coefficient in value.list()]
    )


def interpolate_point_images(samples, numerator_degree, denominator_degree, precision_gate):
    unknown_count = numerator_degree + 1 + denominator_degree
    rows = []
    right_hand_side = []
    used = []
    for point, image in samples:
        rows.append(
            [point**index for index in range(numerator_degree + 1)]
            + [-image * point**index for index in range(denominator_degree)]
        )
        right_hand_side.append(image * point**denominator_degree)
        used.append((point, image))
        if len(rows) == unknown_count:
            break
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
    result = padic_u_field(numerator) / padic_u_field(denominator)
    held_out = 0
    for point, image in samples[unknown_count:]:
        result_denominator = result.denominator()(point)
        if not result_denominator or result_denominator.valuation() != 0:
            continue
        difference = result.numerator()(point) / result_denominator - image
        if difference and difference.valuation() < precision_gate:
            raise ArithmeticError("pointwise Newton interpolation held-out failed")
        held_out += 1
        if held_out == 3:
            break
    if held_out != 3:
        raise ArithmeticError("insufficient pointwise Newton held-outs")
    return result


def pointwise_factor_newton(valuation):
    required = max(
        int(finite_u_field(value).numerator().degree())
        + 1
        + int(finite_u_field(value).denominator().degree())
        for value in [L_finite[0]] + list(Q_finite[:4]) + list(D_finite[:4])
    ) + 3
    samples = []
    for point in hensel_sample_points:
        try:
            local_L = evaluate_w_at(L, point)
            local_Q = evaluate_w_at(Q, point)
            local_D = evaluate_w_at(D, point)
            local_target = evaluate_w_at(target, point)
        except ZeroDivisionError:
            continue
        local_residual = local_L**3 * local_Q**2 * local_D - local_target
        columns = (
            [3 * local_L**2 * local_Q**2 * local_D]
            + [
                2 * local_L**3 * local_Q * local_D * W_checkpoint**index
                for index in range(4)
            ]
            + [
                local_L**3 * local_Q**2 * W_checkpoint**index
                for index in range(4)
            ]
        )
        jacobian = Matrix(
            padic_field,
            9,
            9,
            lambda row, column: columns[column][row],
        )
        determinant = jacobian.det()
        if not determinant or determinant.valuation() != 0:
            continue
        correction = jacobian.solve_right(
            Matrix(padic_field, 9, 1, [-local_residual[index] for index in range(9)])
        ).column(0)
        images = (
            [local_L[0] + correction[0]]
            + [local_Q[index] + correction[1 + index] for index in range(4)]
            + [local_D[index] + correction[5 + index] for index in range(4)]
        )
        samples.append((point, images))
        if len(samples) == required:
            break
    if len(samples) != required:
        raise ArithmeticError("insufficient good points for pointwise factor Newton")
    precision_gate = min(2 * valuation, args.digits)
    finite_values = [L_finite[0]] + list(Q_finite[:4]) + list(D_finite[:4])
    interpolated = []
    for coordinate, finite_value in enumerate(finite_values):
        finite_value = finite_u_field(finite_value)
        interpolated.append(
            interpolate_point_images(
                [(point, images[coordinate]) for point, images in samples],
                int(finite_value.numerator().degree()),
                int(finite_value.denominator().degree()),
                precision_gate,
            )
        )
    return (
        padic_w_ring([interpolated[0], padic_u_field.one()]),
        padic_w_ring(interpolated[1:5] + [padic_u_field.one()]),
        padic_w_ring(interpolated[5:9] + [padic_u_field.one()]),
    )


def pointwise_root_newton(root_residual, valuation):
    finite_values = list(finite_repeated_root)
    required = max(
        int(finite_u_field(value).numerator().degree())
        + 1
        + int(finite_u_field(value).denominator().degree())
        for value in finite_values
    ) + 3
    second_derivative = padic_w_ring(
        padic_cubic.derivative().derivative()(repeated_root)
    )
    samples = []
    for point in hensel_sample_points:
        try:
            local_root = evaluate_w_at(repeated_root, point)
            local_residual = evaluate_w_at(root_residual, point)
            local_second = evaluate_w_at(second_derivative, point)
            local_conductor = evaluate_w_at(conductor, point)
        except ZeroDivisionError:
            continue
        try:
            inverse = local_second.inverse_mod(local_conductor)
        except (ArithmeticError, ZeroDivisionError):
            continue
        corrected = (local_root - local_residual * inverse) % local_conductor
        samples.append((point, list(corrected)))
        if len(samples) == required:
            break
    if len(samples) != required:
        raise ArithmeticError("insufficient good points for pointwise root Newton")
    precision_gate = min(2 * valuation, args.digits)
    interpolated = []
    for coordinate, finite_value in enumerate(finite_values):
        finite_value = finite_u_field(finite_value)
        interpolated.append(
            interpolate_point_images(
                [(point, images[coordinate]) for point, images in samples],
                int(finite_value.numerator().degree()),
                int(finite_value.denominator().degree()),
                precision_gate,
            )
        )
    return padic_w_ring(interpolated)
'''

source = replace_once(
    source,
    "\n\nvaluation_history = []\nwhile True:\n",
    checkpoint_helpers
    + '''

valuation_history = (
    list(checkpoint_payload["factor_valuation_history_prefix"])
    if checkpoint_payload is not None
    else []
)
while True:
''',
    "factor checkpoint helpers",
)
source = replace_once(
    source,
    "    valuation = int(w_polynomial_valuation(residual))\n",
    '''    raw_valuation = w_polynomial_valuation(residual)
    valuation = args.digits if raw_valuation == infinity else int(raw_valuation)
''',
    "infinite terminal factor valuation",
)
source = replace_once(
    source,
    '    print(f"Q80THIRDQ12PADICFACTORS_PROGRESS|stage=factors|valuation={valuation}", flush=True)\n'
    '    if valuation >= args.digits:\n',
    '    print(f"Q80THIRDQ12PADICFACTORS_PROGRESS|stage=factors|valuation={valuation}", flush=True)\n'
    '    write_checkpoint("factors", valuation)\n'
    '    if valuation >= args.digits:\n',
    "factor checkpoint call",
)
source = replace_once(
    source,
    '''    right_hand_side = Matrix(
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
''',
    '''    if args.pointwise_newton:
        L, Q, D = pointwise_factor_newton(valuation)
    else:
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
''',
    "pointwise factor Newton",
)
source = replace_once(
    source,
    '''repeated_root = lift_finite_polynomial(finite_repeated_root)
root_valuation_history = []
while True:
''',
    '''repeated_root = (
    load_checkpoint_w(checkpoint_payload["repeated_root"])
    if checkpoint_payload is not None and checkpoint_payload["stage"] == "root"
    else lift_finite_polynomial(finite_repeated_root)
)
root_valuation_history = (
    list(checkpoint_payload["root_valuation_history_prefix"])
    if checkpoint_payload is not None and checkpoint_payload["stage"] == "root"
    else []
)
while True:
''',
    "root checkpoint load",
)
source = replace_once(
    source,
    "    root_valuation = int(w_polynomial_valuation(root_residual))\n",
    '''    raw_root_valuation = w_polynomial_valuation(root_residual)
    root_valuation = (
        args.digits if raw_root_valuation == infinity else int(raw_root_valuation)
    )
''',
    "infinite terminal root valuation",
)
source = replace_once(
    source,
    '    print(f"Q80THIRDQ12PADICFACTORS_PROGRESS|stage=root|valuation={root_valuation}", flush=True)\n'
    '    if root_valuation >= args.digits:\n',
    '    print(f"Q80THIRDQ12PADICFACTORS_PROGRESS|stage=root|valuation={root_valuation}", flush=True)\n'
    '    write_checkpoint("root", root_valuation)\n'
    '    if root_valuation >= args.digits:\n',
    "root checkpoint call",
)
source = replace_once(
    source,
    '''    finite_residual = finite_w_ring(
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
''',
    '''    if args.pointwise_newton:
        repeated_root = pointwise_root_newton(root_residual, root_valuation)
    else:
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
''',
    "pointwise root Newton",
)
source = replace_once(
    source,
    '''compressed_factor_residual_valuation = int(
    w_polynomial_valuation(L**3 * Q**2 * D - target)
)
''',
    '''raw_compressed_factor_residual_valuation = w_polynomial_valuation(
    L**3 * Q**2 * D - target
)
compressed_factor_residual_valuation = (
    args.digits
    if raw_compressed_factor_residual_valuation == infinity
    else int(raw_compressed_factor_residual_valuation)
)
''',
    "infinite compressed factor valuation",
)

namespace = {
    "__file__": str(CANONICAL),
    "__name__": "__main__",
    "__checkpoint_runner_path__": str(Path(__file__).resolve()),
}
exec(compile(source, str(CANONICAL), "exec"), namespace)
