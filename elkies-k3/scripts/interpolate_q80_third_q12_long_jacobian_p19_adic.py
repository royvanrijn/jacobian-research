#!/usr/bin/env python3
"""Interpolate the exact-gauge long Jacobian over (Z/19^5)[omega](U)."""

import argparse
import glob
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_TRANSPORT = RESULTS / "q80-third-q12-long-jacobians-exact-quadratic-gauge.json"
DEFAULT_OUTPUT = RESULTS / "q80-third-q12-long-jacobian-p19-adic-precision5.json"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--transport", type=Path, default=DEFAULT_TRANSPORT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
args.transport = args.transport.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


prime = 19
digits = 5
modulus = prime**digits
source = json.loads(
    (RESULTS / "q80-third-q12-exact-pencil-p19-adic-precision64.json").read_text()
)
omega_square = int(source["quadratic_field"]["omega_square_modulus"]) % modulus
ZERO = (0, 0)
ONE = (1, 0)


def c(value):
    return int(value[0]) % modulus, int(value[1]) % modulus


def add(left, right):
    return (left[0] + right[0]) % modulus, (left[1] + right[1]) % modulus


def neg(value):
    return -value[0] % modulus, -value[1] % modulus


def sub(left, right):
    return add(left, neg(right))


def mul(left, right):
    return (
        (left[0] * right[0] + omega_square * left[1] * right[1]) % modulus,
        (left[0] * right[1] + left[1] * right[0]) % modulus,
    )


def power(value, exponent):
    result = ONE
    base = value
    while exponent:
        if exponent & 1:
            result = mul(result, base)
        base = mul(base, base)
        exponent >>= 1
    return result


def is_unit(value):
    return (value[0] * value[0] - omega_square * value[1] * value[1]) % prime != 0


def inverse(value):
    norm = (value[0] * value[0] - omega_square * value[1] * value[1]) % modulus
    if norm % prime == 0:
        raise ZeroDivisionError("non-unit quadratic coefficient")
    inverse_norm = pow(norm, -1, modulus)
    return value[0] * inverse_norm % modulus, -value[1] * inverse_norm % modulus


def divide(left, right):
    return mul(left, inverse(right))


sample_paths = sorted(
    Path(path).resolve()
    for path in glob.glob(str(RESULTS / "q80-third-q12-p19-adic-U*.json"))
)
samples_by_base = {}
for path in sample_paths:
    payload = json.loads(path.read_text())
    if payload.get("status") != "PASS_EXACT_THIRD_Q12_WEIERSTRASS_P19_ADIC_SAMPLE":
        continue
    base = tuple(payload["specialization"]["base_U_coefficients_1_omega"])
    if base in samples_by_base:
        raise ArithmeticError(f"duplicate p-adic sample at U={base}")
    samples_by_base[base] = (path, payload)
samples = [samples_by_base[key] for key in sorted(samples_by_base, key=lambda pair: (pair[1], pair[0]))]
if len(samples) < 19:
    raise ArithmeticError(f"need at least 19 p-adic samples, found {len(samples)}")
training = samples[:-2]
held_out = samples[-2:]
if len(training) < 17:
    raise ArithmeticError("need seventeen training samples for the degree 8/8 coefficient")


def solve_overdetermined(rows, unknown_count):
    matrix = [[c(value) for value in row] for row in rows]
    rank = 0
    pivot_rows = []
    for column in range(unknown_count):
        pivot = next(
            (index for index in range(rank, len(matrix)) if is_unit(matrix[index][column])),
            None,
        )
        if pivot is None:
            raise ArithmeticError(f"interpolation matrix lost a unit pivot at column {column}")
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = inverse(matrix[rank][column])
        matrix[rank] = [mul(value, scale) for value in matrix[rank]]
        for row_index in range(len(matrix)):
            if row_index == rank:
                continue
            factor = matrix[row_index][column]
            if factor == ZERO:
                continue
            matrix[row_index] = [
                sub(matrix[row_index][index], mul(factor, matrix[rank][index]))
                for index in range(unknown_count + 1)
            ]
        pivot_rows.append(rank)
        rank += 1
    for row in matrix[rank:]:
        if any(value != ZERO for value in row[:unknown_count]) or row[-1] != ZERO:
            raise ArithmeticError("overdetermined interpolation equations are inconsistent")
    return [matrix[index][-1] for index in pivot_rows]


names = ("a1", "a2", "a3", "a4", "a6")
degree_bounds = dict(zip(names, (2, 4, 4, 6, 8)))
interpolated = {}
for coefficient_index, name in enumerate(names):
    degree = degree_bounds[name]
    unknown_count = 2 * degree + 1
    rows = []
    for unused_path, payload in training:
        u_value = c(payload["specialization"]["base_U_coefficients_1_omega"])
        y_value = c(payload["weierstrass"]["a1_a2_a3_a4_a6_mod_19_power_1_omega"][coefficient_index])
        powers = [power(u_value, exponent) for exponent in range(degree + 1)]
        rows.append(
            powers
            + [neg(mul(y_value, powers[exponent])) for exponent in range(degree)]
            + [mul(y_value, powers[degree])]
        )
    solution = solve_overdetermined(rows, unknown_count)
    numerator = solution[: degree + 1]
    denominator = solution[degree + 1 :] + [ONE]
    held_out_checks = []
    for path, payload in held_out:
        u_value = c(payload["specialization"]["base_U_coefficients_1_omega"])
        expected = c(payload["weierstrass"]["a1_a2_a3_a4_a6_mod_19_power_1_omega"][coefficient_index])
        numerator_value = ZERO
        denominator_value = ZERO
        for exponent, value in enumerate(numerator):
            numerator_value = add(numerator_value, mul(value, power(u_value, exponent)))
        for exponent, value in enumerate(denominator):
            denominator_value = add(denominator_value, mul(value, power(u_value, exponent)))
        if not is_unit(denominator_value):
            raise ArithmeticError(f"held-out denominator for {name} is not a unit at {u_value}")
        actual = divide(numerator_value, denominator_value)
        if actual != expected:
            raise ArithmeticError(f"held-out p-adic replay failed for {name} at {u_value}")
        held_out_checks.append({"path": str(path.relative_to(ROOT)), "base_U": list(u_value)})
    interpolated[name] = {
        "degrees_numerator_denominator": [degree, degree],
        "numerator_coefficients_low_to_high_U_mod_19_power_1_omega": [list(value) for value in numerator],
        "denominator_coefficients_low_to_high_U_mod_19_power_1_omega": [list(value) for value in denominator],
        "training_equation_count": len(training),
        "held_out_replays": held_out_checks,
    }

transport = json.loads(args.transport.read_text())
transport_model = transport["transported_models"]["19"]["weierstrass"]
for name in names:
    p_adic = interpolated[name]
    finite_record = transport_model[name]
    for p_adic_key, finite_key in (
        ("numerator_coefficients_low_to_high_U_mod_19_power_1_omega", "numerator_coefficients_low_to_high_1_omega"),
        ("denominator_coefficients_low_to_high_U_mod_19_power_1_omega", "denominator_coefficients_low_to_high_1_omega"),
    ):
        reduction = [[value[0] % prime, value[1] % prime] for value in p_adic[p_adic_key]]
        if reduction != finite_record[finite_key]:
            raise ArithmeticError(f"interpolated {name} does not reduce to transported p=19 model")


def interpolate_scalar_function(value_getter, finite_record, label):
    numerator_degree, denominator_degree = finite_record["degrees_numerator_denominator"]
    if numerator_degree == -1:
        for unused_path, payload in map_samples:
            if c(value_getter(payload)) != ZERO:
                raise ArithmeticError(f"structural zero map slot became nonzero: {label}")
        return {
            "degrees_numerator_denominator": [-1, 0],
            "numerator_coefficients_low_to_high_U_mod_19_power_1_omega": [],
            "denominator_coefficients_low_to_high_U_mod_19_power_1_omega": [[1, 0]],
            "training_equation_count": len(map_training),
            "held_out_replay_count": len(map_held_out),
        }
    if numerator_degree != denominator_degree:
        raise ArithmeticError(f"non-square rational degree shape for {label}")
    degree = numerator_degree
    unknown_count = 2 * degree + 1
    rows = []
    for unused_path, payload in map_training:
        u_value = c(payload["specialization"]["base_U_coefficients_1_omega"])
        y_value = c(value_getter(payload))
        powers = [power(u_value, exponent) for exponent in range(degree + 1)]
        rows.append(
            powers
            + [neg(mul(y_value, powers[exponent])) for exponent in range(degree)]
            + [mul(y_value, powers[degree])]
        )
    try:
        solution = solve_overdetermined(rows, unknown_count)
    except ArithmeticError as error:
        raise ArithmeticError(f"{label}: {error}") from error
    numerator = solution[: degree + 1]
    denominator = solution[degree + 1 :] + [ONE]
    for unused_path, payload in map_held_out:
        u_value = c(payload["specialization"]["base_U_coefficients_1_omega"])
        expected = c(value_getter(payload))
        numerator_value = ZERO
        denominator_value = ZERO
        u_power = ONE
        for exponent in range(degree + 1):
            numerator_value = add(numerator_value, mul(numerator[exponent], u_power))
            denominator_value = add(denominator_value, mul(denominator[exponent], u_power))
            u_power = mul(u_power, u_value)
        if not is_unit(denominator_value) or divide(numerator_value, denominator_value) != expected:
            raise ArithmeticError(f"held-out map replay failed for {label}")
    for values, finite_key in (
        (numerator, "numerator_coefficients_low_to_high_1_omega"),
        (denominator, "denominator_coefficients_low_to_high_1_omega"),
    ):
        reduction = [[value[0] % prime, value[1] % prime] for value in values]
        if reduction != finite_record[finite_key]:
            raise ArithmeticError(f"interpolated map slot does not reduce to p=19: {label}")
    return {
        "degrees_numerator_denominator": [degree, degree],
        "numerator_coefficients_low_to_high_U_mod_19_power_1_omega": [list(value) for value in numerator],
        "denominator_coefficients_low_to_high_U_mod_19_power_1_omega": [list(value) for value in denominator],
        "training_equation_count": len(map_training),
        "held_out_replay_count": len(map_held_out),
    }


map_transport_path = RESULTS / "q80-third-q12-birational-maps-exact-quadratic-gauge.json"
map_transport = json.loads(map_transport_path.read_text())
finite_maps = map_transport["transported_maps"]["19"]


def reduce_sample_rational_W(numerator, denominator, numerator_degree, denominator_degree, label):
    numerator = [c(value) for value in numerator]
    denominator = [c(value) for value in denominator]
    unknown_count = numerator_degree + 1 + denominator_degree
    maximum_equation_degree = max(
        len(numerator) - 1 + denominator_degree,
        len(denominator) - 1 + numerator_degree,
    )
    rows = []
    for total_degree in range(maximum_equation_degree + 1):
        row = []
        # Coefficients of reduced N enter through -denominator*N.
        for degree in range(numerator_degree + 1):
            index = total_degree - degree
            row.append(neg(denominator[index]) if 0 <= index < len(denominator) else ZERO)
        # Non-leading coefficients of reduced D enter through numerator*D.
        for degree in range(denominator_degree):
            index = total_degree - degree
            row.append(numerator[index] if 0 <= index < len(numerator) else ZERO)
        leading_index = total_degree - denominator_degree
        leading_term = numerator[leading_index] if 0 <= leading_index < len(numerator) else ZERO
        rows.append(row + [neg(leading_term)])
    try:
        solution = solve_overdetermined(rows, unknown_count)
    except ArithmeticError as error:
        raise ArithmeticError(f"W-reduction failed for {label}: {error}") from error
    return solution[: numerator_degree + 1], solution[numerator_degree + 1 :] + [ONE]


normalized_forward_samples = {}
map_samples = []
rejected_map_samples = []
for path, payload in samples:
    base = tuple(payload["specialization"]["base_U_coefficients_1_omega"])
    candidate = {}
    try:
        for coordinate, sample_key in (
            ("X", "X_weierstrass_as_function_of_W_z"),
            ("Y", "Y_weierstrass_as_function_of_W_z"),
        ):
            candidate[coordinate] = []
            sample_components = payload["birational_maps"]["forward"][sample_key][
                "coefficients_low_to_high_z"
            ]
            for component_index, (sample_component, finite_component) in enumerate(
                zip(sample_components, finite_maps["forward_long"][coordinate])
            ):
                numerator_degree, denominator_degree = finite_component[
                    "degrees_W_numerator_denominator"
                ]
                reduced = reduce_sample_rational_W(
                    sample_component["numerator_coefficients_low_to_high_W_mod_19_power_1_omega"],
                    sample_component["denominator_coefficients_low_to_high_W_mod_19_power_1_omega"],
                    numerator_degree,
                    denominator_degree,
                    f"U={base}.{coordinate}.z{component_index}",
                )
                candidate[coordinate].append(reduced)
    except ArithmeticError as error:
        rejected_map_samples.append(
            {"path": str(path.relative_to(ROOT)), "base_U": list(base), "reason": str(error)}
        )
        continue
    normalized_forward_samples[base] = candidate
    map_samples.append((path, payload))
if len(map_samples) < 83:
    raise ArithmeticError(
        f"need 83 canonically reduced map samples, found {len(map_samples)}"
    )
map_training = map_samples[:-2]
map_held_out = map_samples[-2:]

interpolated_forward = {}
for coordinate, sample_key in (
    ("X", "X_weierstrass_as_function_of_W_z"),
    ("Y", "Y_weierstrass_as_function_of_W_z"),
):
    interpolated_forward[coordinate] = []
    finite_components = finite_maps["forward_long"][coordinate]
    for component_index, finite_component in enumerate(finite_components):
        component_output = {"coefficients_low_to_high_z_component": component_index}
        for kind in ("numerator", "denominator"):
            finite_records = finite_component[f"{kind}_coefficients_low_to_high_W"]
            sample_record_key = f"{kind}_coefficients_low_to_high_W_mod_19_power_1_omega"
            records = []
            for w_index, finite_record in enumerate(finite_records):
                def getter(payload, component_index=component_index, w_index=w_index):
                    base = tuple(payload["specialization"]["base_U_coefficients_1_omega"])
                    reduced = normalized_forward_samples[base][coordinate][component_index]
                    values = reduced[0 if kind == "numerator" else 1]
                    return values[w_index] if w_index < len(values) else [0, 0]

                records.append(
                    interpolate_scalar_function(
                        getter,
                        finite_record,
                        f"forward.{coordinate}.z{component_index}.{kind}.W{w_index}",
                    )
                )
            component_output[f"{kind}_coefficients_low_to_high_W"] = records
        interpolated_forward[coordinate].append(component_output)

interpolated_inverse = {}
for transport_target, sample_target in (("W", "W"), ("old_x", "z")):
    finite_target = finite_maps["inverse_long"][transport_target]
    target_output = {
        "weighted_bound": finite_target["weighted_bound"],
        "monomials_X_power_Y_power": finite_target["monomials_X_power_Y_power"],
        "formula": finite_target["formula"],
    }
    for kind in ("numerator", "denominator"):
        finite_records = finite_target[f"{kind}_coefficients"]
        sample_key = f"{kind}_coefficients_mod_19_power_1_omega"
        records = []
        for monomial_index, finite_record in enumerate(finite_records):
            def getter(payload, monomial_index=monomial_index):
                return payload["birational_maps"]["inverse"][sample_target][sample_key][monomial_index]

            records.append(
                interpolate_scalar_function(
                    getter,
                    finite_record,
                    f"inverse.{transport_target}.{kind}.m{monomial_index}",
                )
            )
        target_output[f"{kind}_coefficients"] = records
    interpolated_inverse[transport_target] = target_output

interpolated_maps = {
    "forward_long": interpolated_forward,
    "inverse_long": interpolated_inverse,
    "all_scalar_slots_reduce_to_transported_p19": True,
    "all_scalar_slots_replay_two_held_out_padic_samples": True,
    "sample_counts": {
        "accepted": len(map_samples),
        "training": len(map_training),
        "held_out": len(map_held_out),
        "rejected_during_canonical_W_reduction": len(rejected_map_samples),
    },
    "rejected_samples": rejected_map_samples,
}

output = {
    "schema": "elkies-k3.q80-third-q12-long-jacobian-p19-adic.v1",
    "status": "PASS_EXACT_THIRD_Q12_LONG_CHILD_INTERPOLATION_P19_ADIC",
    "specialization": {"u": "-2", "prime": prime, "digits": digits, "modulus": modulus},
    "weierstrass": interpolated,
    "birational_maps": interpolated_maps,
    "samples": {
        "training": [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "base_U": payload["specialization"]["base_U_coefficients_1_omega"],
            }
            for path, payload in training
        ],
        "held_out": [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "base_U": payload["specialization"]["base_U_coefficients_1_omega"],
            }
            for path, payload in held_out
        ],
    },
    "positive_control": {
        "path": str(args.transport.relative_to(ROOT)),
        "sha256": sha256(args.transport),
        "literal_reduction_mod_19_of_all_coefficients": True,
        "map_transport": {
            "path": str(map_transport_path.relative_to(ROOT)),
            "sha256": sha256(map_transport_path),
            "literal_reduction_mod_19_of_all_map_slots": True,
        },
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "the complete exact-gauge long Weierstrass equation modulo 19^5 as rational functions of U",
            "two held-out p-adic sample replays for every long coefficient",
            "literal reduction to the independently transported generic p=19 long model",
            "the complete forward and inverse maps modulo 19^5 as rational functions of U",
            "two held-out p-adic replays and transported mod-19 replay for every map scalar function",
        ],
        "not_proved": [
            "characteristic-zero reconstruction, minimization, fibres, or marking",
        ],
    },
    "reproduce": "python3 elkies-k3/scripts/interpolate_q80_third_q12_long_jacobian_p19_adic.py",
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"p-adic long-Jacobian interpolation artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    "Q80THIRDQ12PADICLONG|p=19|digits=5|samples={}+{}|degrees=2,4,4,6,8|maps=both|"
    "heldout=PASS|transport=PASS|status=PASS_EXACT_THIRD_Q12_LONG_CHILD_INTERPOLATION_P19_ADIC".format(
        len(training), len(held_out)
    )
)
