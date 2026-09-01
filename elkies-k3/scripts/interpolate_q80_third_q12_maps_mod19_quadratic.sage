#!/usr/bin/env sage -python
"""Interpolate and certify the generic p=19 third-q12 birational maps.

Solve jointly for compact bivariate numerators and denominators, replay
held-out fibres, and verify both directions in the generic function field.
"""

import hashlib
import json
from pathlib import Path

from sage.all import FunctionField, GF, Matrix, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
BATCH = ROOT / "artifacts/generated-results/q80-third-q12-p19-weierstrass-sample-batch.json"
PENCIL = ROOT / "artifacts/generated-results/q80-third-q12-um2-p19-resolved-pencil.json"
JACOBIAN = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-interpolated.json"
MINIMAL = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-minimal.json"
OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-birational-maps.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


batch = json.loads(BATCH.read_text())
pencil = json.loads(PENCIL.read_text())
jacobian = json.loads(JACOBIAN.read_text())
minimal = json.loads(MINIMAL.read_text())
expected_inputs = (
    (batch, "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_MOD19_QUADRATIC"),
    (pencil, "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC"),
    (jacobian, "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_MOD19_QUADRATIC"),
    (minimal, "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_MOD19_QUADRATIC"),
)
if any(payload.get("status") != status for payload, status in expected_inputs):
    raise ValueError("one or more map inputs are not certified")

base_finite = GF(19)
modulus_ring = PolynomialRing(base_finite, "m")
m = modulus_ring.gen()
finite = GF(19**2, "r", modulus=m**2 + 12 * m + 3)
r = finite.gen()
v_ring = PolynomialRing(finite, "V")
V = v_ring.gen()
v_field = v_ring.fraction_field()
w_test_ring = PolynomialRing(finite, "Wtest")


def element(coordinates):
    return finite(coordinates[0]) + finite(coordinates[1]) * r


def load_sample(summary):
    path = ROOT / summary["path"]
    if sha256(path) != summary["sha256"]:
        raise ArithmeticError(f"sample hash mismatch: {summary['path']}")
    payload = json.loads(path.read_text())
    if payload.get("status") != "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_MOD19_QUADRATIC":
        raise ArithmeticError(f"uncertified sample: {summary['path']}")
    return element(summary["new_base_coefficients_1_r"]), payload


training = [load_sample(summary) for summary in batch["training_samples"]]
held_out = [load_sample(summary) for summary in batch["held_out_samples"]]
all_samples = training + held_out
if len(training) < 100 or len(held_out) != 8:
    raise ArithmeticError("insufficient training/held-out split")
solve_samples = training[:32]
replay_samples = training[32:] + held_out


def forward_local(sample, coordinate, old_x_power):
    record = sample["birational_maps"]["forward"][
        f"{coordinate}_weierstrass_as_function_of_W_old_x"
    ]["coefficients_low_to_high_old_x"][old_x_power]
    return (
        [element(value) for value in record["numerator_coefficients_low_to_high_W"]],
        [element(value) for value in record["denominator_coefficients_low_to_high_W"]],
    )


def evaluate_polynomial_list(polynomials, value):
    return w_test_ring([polynomial(value) for polynomial in polynomials])


def solve_forward_block(coordinate, old_x_power, numerator_w_degree, denominator_w_degree):
    for total_v_degree in range(25):
        for numerator_v_degree in range(total_v_degree + 1):
            denominator_v_degree = total_v_degree - numerator_v_degree
            column_count = (numerator_w_degree + 1) * (numerator_v_degree + 1)
            column_count += (denominator_w_degree + 1) * (denominator_v_degree + 1)
            rows = []
            for value, sample in solve_samples:
                local_numerator, local_denominator = forward_local(sample, coordinate, old_x_power)
                coefficient_count = max(
                    numerator_w_degree + len(local_denominator),
                    denominator_w_degree + len(local_numerator),
                )
                for w_power in range(coefficient_count):
                    row = []
                    for global_w_power in range(numerator_w_degree + 1):
                        index = w_power - global_w_power
                        coefficient = local_denominator[index] if 0 <= index < len(local_denominator) else finite.zero()
                        row.extend(coefficient * value**i for i in range(numerator_v_degree + 1))
                    for global_w_power in range(denominator_w_degree + 1):
                        index = w_power - global_w_power
                        coefficient = local_numerator[index] if 0 <= index < len(local_numerator) else finite.zero()
                        row.extend(-coefficient * value**i for i in range(denominator_v_degree + 1))
                    if any(row):
                        rows.append(row)
            kernel = Matrix(finite, rows, ncols=column_count).right_kernel()
            if kernel.dimension() != 1:
                continue
            vector = kernel.basis()[0]
            vector /= next(coefficient for coefficient in vector if coefficient)
            cursor = 0
            numerator = []
            for unused in range(numerator_w_degree + 1):
                numerator.append(v_ring(list(vector[cursor : cursor + numerator_v_degree + 1])))
                cursor += numerator_v_degree + 1
            denominator = []
            for unused in range(denominator_w_degree + 1):
                denominator.append(v_ring(list(vector[cursor : cursor + denominator_v_degree + 1])))
                cursor += denominator_v_degree + 1
            for value, sample in replay_samples:
                local_numerator, local_denominator = forward_local(sample, coordinate, old_x_power)
                global_numerator = evaluate_polynomial_list(numerator, value)
                global_denominator = evaluate_polynomial_list(denominator, value)
                if (
                    not global_denominator
                    or global_numerator * w_test_ring(local_denominator)
                    != global_denominator * w_test_ring(local_numerator)
                ):
                    break
            else:
                return {
                    "numerator": numerator,
                    "denominator": denominator,
                    "degrees_V_numerator_denominator": [numerator_v_degree, denominator_v_degree],
                    "degrees_W_numerator_denominator": [numerator_w_degree, denominator_w_degree],
                }
    raise ArithmeticError(f"no joint forward block for {coordinate}, old_x^{old_x_power}")


forward_shapes = {
    "X": ((6, 4), (4, 5), (1, 5)),
    "Y": ((7, 4), (5, 5), (2, 5)),
}
forward_records = {
    coordinate: [
        solve_forward_block(coordinate, old_x_power, *shape)
        for old_x_power, shape in enumerate(shapes)
    ]
    for coordinate, shapes in forward_shapes.items()
}

inverse_labels = {
    "W": [(0, 0), (1, 0), (2, 0), (0, 1)],
    "old_x": [(power, 0) for power in range(6)] + [(power, 1) for power in range(4)],
}


def inverse_local(sample, target):
    record = sample["birational_maps"]["inverse"][target]
    local_labels = [tuple(label) for label in record["monomials_X_power_Y_power"]]
    numerator = dict(zip(local_labels, map(element, record["numerator_coefficients"])))
    denominator = dict(zip(local_labels, map(element, record["denominator_coefficients"])))
    labels = inverse_labels[target]
    return (
        [numerator.get(label, finite.zero()) for label in labels],
        [denominator.get(label, finite.zero()) for label in labels],
    )


def sample_weierstrass_coefficients(sample):
    return [element(value) for value in sample["weierstrass"]["a1_a2_a3_a4_a6"]]


def multiply_monomials(left, right, coefficients):
    x_power = left[0] + right[0]
    y_power = left[1] + right[1]
    if y_power < 2:
        return {(x_power, y_power): finite.one()}
    a1, a2, a3, a4, a6 = coefficients
    return {
        (x_power + 3, 0): finite.one(),
        (x_power + 2, 0): a2,
        (x_power + 1, 0): a4,
        (x_power, 0): a6,
        (x_power + 1, 1): -a1,
        (x_power, 1): -a3,
    }


def multiply_functions(left, right, coefficients):
    result = {}
    for left_label, left_value in left.items():
        for right_label, right_value in right.items():
            for label, structure_constant in multiply_monomials(left_label, right_label, coefficients).items():
                result[label] = result.get(label, finite.zero()) + left_value * right_value * structure_constant
    return {label: value for label, value in result.items() if value}


def solve_inverse_target(target):
    labels = inverse_labels[target]
    for total_v_degree in range(25):
        for numerator_v_degree in range(total_v_degree + 1):
            denominator_v_degree = total_v_degree - numerator_v_degree
            numerator_columns = len(labels) * (numerator_v_degree + 1)
            column_count = numerator_columns + len(labels) * (denominator_v_degree + 1)
            rows = []
            for value, sample in solve_samples:
                local_numerator, local_denominator = inverse_local(sample, target)
                coefficients = sample_weierstrass_coefficients(sample)
                row_by_output = {}
                for global_index, global_label in enumerate(labels):
                    for local_index, local_label in enumerate(labels):
                        for is_numerator, local_value, v_degree, offset, sign in (
                            (True, local_denominator[local_index], numerator_v_degree, 0, 1),
                            (False, local_numerator[local_index], denominator_v_degree, numerator_columns, -1),
                        ):
                            if not local_value:
                                continue
                            for output_label, structure_constant in multiply_monomials(
                                global_label, local_label, coefficients
                            ).items():
                                row = row_by_output.setdefault(output_label, [finite.zero()] * column_count)
                                for v_power in range(v_degree + 1):
                                    column = offset + global_index * (v_degree + 1) + v_power
                                    row[column] += sign * local_value * structure_constant * value**v_power
                rows.extend(row for row in row_by_output.values() if any(row))
            kernel = Matrix(finite, rows, ncols=column_count).right_kernel()
            if kernel.dimension() != 1:
                continue
            vector = kernel.basis()[0]
            vector /= next(coefficient for coefficient in vector if coefficient)
            cursor = 0
            numerator = []
            for unused in labels:
                numerator.append(v_ring(list(vector[cursor : cursor + numerator_v_degree + 1])))
                cursor += numerator_v_degree + 1
            denominator = []
            for unused in labels:
                denominator.append(v_ring(list(vector[cursor : cursor + denominator_v_degree + 1])))
                cursor += denominator_v_degree + 1
            for value, sample in replay_samples:
                local_numerator, local_denominator = inverse_local(sample, target)
                coefficients = sample_weierstrass_coefficients(sample)
                global_numerator = {
                    label: polynomial(value) for label, polynomial in zip(labels, numerator) if polynomial(value)
                }
                global_denominator = {
                    label: polynomial(value) for label, polynomial in zip(labels, denominator) if polynomial(value)
                }
                left = multiply_functions(global_numerator, dict(zip(labels, local_denominator)), coefficients)
                right = multiply_functions(global_denominator, dict(zip(labels, local_numerator)), coefficients)
                if left != right or not global_denominator:
                    break
            else:
                return {
                    "numerator": numerator,
                    "denominator": denominator,
                    "degrees_V_numerator_denominator": [numerator_v_degree, denominator_v_degree],
                }
    raise ArithmeticError(f"no joint inverse for {target}")


# Literal generic replay.
old_function = FunctionField(v_field, "W")
W = old_function.gen()
x_ring = PolynomialRing(old_function, "zpoly")
zpoly = x_ring.gen()
resolved_equation = x_ring.zero()
for t_degree, w_degree, x_degree, coordinates in pencil["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    resolved_equation += v_field(element(coordinates) * V**t_degree) * W**w_degree * zpoly**x_degree
if resolved_equation.degree() != 3:
    raise ArithmeticError("generic resolved equation is not cubic")
# Irreducibility is certified by any retained degree-preserving specialization:
# a generic factorization would specialize away from finitely many poles, while
# every retained worker sample proves its specialized cubic irreducible.
curve_function = old_function.extension(resolved_equation.monic(), "old_x")
old_x = curve_function.gen()


def bivariate_function(record):
    numerator = sum(
        (old_function(v_field(polynomial)) * W**power for power, polynomial in enumerate(record["numerator"])),
        old_function.zero(),
    )
    denominator = sum(
        (old_function(v_field(polynomial)) * W**power for power, polynomial in enumerate(record["denominator"])),
        old_function.zero(),
    )
    return curve_function(numerator / denominator)


X_long = sum(
    (bivariate_function(record) * old_x**power for power, record in enumerate(forward_records["X"])),
    curve_function.zero(),
)
Y_long = sum(
    (bivariate_function(record) * old_x**power for power, record in enumerate(forward_records["Y"])),
    curve_function.zero(),
)


def polynomial_from_record(coordinates):
    return v_ring([element(value) for value in coordinates])


def rational_from_record(record):
    return v_field(
        polynomial_from_record(record["numerator_coefficients_low_to_high_1_r"])
        / polynomial_from_record(record["denominator_coefficients_low_to_high_1_r"])
    )


a1, a2, a3, a4, a6 = [
    curve_function(rational_from_record(jacobian["weierstrass"][name]))
    for name in ("a1", "a2", "a3", "a4", "a6")
]
if Y_long**2 + a1 * X_long * Y_long + a3 * Y_long != X_long**3 + a2 * X_long**2 + a4 * X_long + a6:
    raise ArithmeticError("generic forward map misses the long Weierstrass equation")


def generic_constant_relation_kernel(function_values):
    rows = []
    for component in range(3):
        values = []
        for value in function_values:
            coefficients = list(curve_function(value).list())
            coefficients += [old_function.zero()] * (3 - len(coefficients))
            values.append(coefficients[component])
        polynomial_ring = values[0].denominator().parent()
        common_denominator = polynomial_ring.one()
        for value in values:
            common_denominator = common_denominator.lcm(value.denominator())
        numerators = []
        for value in values:
            cleared = value * common_denominator
            if cleared.denominator() != 1:
                raise ArithmeticError("failed to clear generic inverse denominator")
            numerators.append(cleared.numerator())
        maximum_degree = max(numerator.degree() for numerator in numerators)
        for degree in range(maximum_degree + 1):
            rows.append([numerator[degree] for numerator in numerators])
    return Matrix(v_field, rows).right_kernel()


def generic_inverse_formula(target, minimum_bound, maximum_bound):
    for bound in range(minimum_bound, maximum_bound + 1):
        labels = []
        monomials = []
        for y_power in range(2):
            for x_power in range(bound // 2 + 1):
                if 2 * x_power + 3 * y_power > bound:
                    continue
                labels.append((x_power, y_power))
                monomials.append(X_long**x_power * Y_long**y_power)
        kernel = generic_constant_relation_kernel(
            tuple(monomials) + tuple(target * monomial for monomial in monomials)
        )
        for relation in kernel.basis():
            split = len(monomials)
            numerator = sum(
                (relation[index] * monomial for index, monomial in enumerate(monomials)),
                curve_function.zero(),
            )
            denominator = sum(
                (relation[split + index] * monomial for index, monomial in enumerate(monomials)),
                curve_function.zero(),
            )
            if denominator and target * denominator + numerator == 0:
                scale = next(coefficient for coefficient in relation if coefficient)
                relation /= scale
                return {
                    "weighted_bound": bound,
                    "labels": labels,
                    "numerator": list(relation[:split]),
                    "denominator": list(relation[split:]),
                }
    raise ArithmeticError(f"no generic inverse through weighted bound {maximum_bound}")


inverse_records = {
    "W": generic_inverse_formula(curve_function(W), 4, 8),
    "old_x": generic_inverse_formula(old_x, 10, 12),
}


def inverse_expression(target):
    labels = inverse_records[target]["labels"]
    record = inverse_records[target]
    numerator = sum(
        (curve_function(coefficient) * X_long**label[0] * Y_long**label[1] for label, coefficient in zip(labels, record["numerator"])),
        curve_function.zero(),
    )
    denominator = sum(
        (curve_function(coefficient) * X_long**label[0] * Y_long**label[1] for label, coefficient in zip(labels, record["denominator"])),
        curve_function.zero(),
    )
    return numerator, denominator


W_numerator, W_denominator = inverse_expression("W")
if curve_function(W) * W_denominator + W_numerator != 0:
    raise ArithmeticError("generic inverse formula for W failed")
old_x_numerator, old_x_denominator = inverse_expression("old_x")
if old_x * old_x_denominator + old_x_numerator != 0:
    raise ArithmeticError("generic inverse formula for old_x failed")


def coordinates(value):
    values = list(finite(value).list()) + [base_finite.zero(), base_finite.zero()]
    return [int(values[0]), int(values[1])]


def polynomial_record(poly):
    return [coordinates(value) for value in poly.list()]


def rational_record(value):
    value = v_field(value)
    return {
        "numerator_coefficients_low_to_high_1_r": polynomial_record(value.numerator()),
        "denominator_coefficients_low_to_high_1_r": polynomial_record(value.denominator()),
        "degrees_numerator_denominator": [int(value.numerator().degree()), int(value.denominator().degree())],
    }


def serialize_joint(record):
    return {
        "numerator_coefficients_low_to_high_auxiliary_power_then_V": [
            polynomial_record(poly) for poly in record["numerator"]
        ],
        "denominator_coefficients_low_to_high_auxiliary_power_then_V": [
            polynomial_record(poly) for poly in record["denominator"]
        ],
        **{key: value for key, value in record.items() if key.startswith("degrees_")},
    }


output = {
    "schema": "elkies-k3.q80-third-q12-birational-maps-modp2.v1",
    "status": "PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_MOD19_QUADRATIC",
    "specialization": {"u": "-2", "prime": 19, "extension_modulus": "r^2+12*r+3"},
    "source_curve": {
        "equation": "sum c[t,w,i]*V^t*W^w*old_x^i=0",
        "coefficient_source": str(PENCIL.relative_to(ROOT)),
    },
    "target_long_weierstrass": {
        "equation": "Y^2+a1(V)XY+a3(V)Y=X^3+a2(V)X^2+a4(V)X+a6(V)",
        "coefficient_source": str(JACOBIAN.relative_to(ROOT)),
    },
    "forward_long": {
        "formula": "X or Y=sum_i N_i(V,W)/D_i(V,W)*old_x^i",
        "old_x_power_support": [0, 1, 2],
        "X": [serialize_joint(record) for record in forward_records["X"]],
        "Y": [serialize_joint(record) for record in forward_records["Y"]],
    },
    "inverse_long": {
        target: {
            "formula": "target=-numerator/denominator",
            "weighted_bound": record["weighted_bound"],
            "monomials_X_power_Y_power": [list(label) for label in record["labels"]],
            "numerator_coefficients": [rational_record(value) for value in record["numerator"]],
            "denominator_coefficients": [rational_record(value) for value in record["denominator"]],
        }
        for target, record in inverse_records.items()
    },
    "long_to_minimal": minimal["long_to_minimal_map"],
    "forward_minimal": {
        "X_min": "q(V)^2*(X_long+b2(V)/12)",
        "Y_min": "q(V)^3*(Y_long+(a1(V)*X_long+a3(V))/2)",
    },
    "inverse_minimal": {
        "X_long": "X_min/q(V)^2-b2(V)/12",
        "Y_long": "Y_min/q(V)^3-(a1(V)*X_long+a3(V))/2",
        "then": "apply inverse_long formulas for W and old_x",
    },
    "validation": {
        "training_samples": len(training),
        "kernel_solve_samples": len(solve_samples),
        "post_solve_replay_samples": len(replay_samples),
        "held_out_samples": len(held_out),
        "joint_cross_multiplication_on_all_samples": len(all_samples),
        "generic_resolved_cubic_irreducible": True,
        "generic_forward_weierstrass_identity": True,
        "generic_inverse_W_identity": True,
        "generic_inverse_old_x_identity": True,
    },
    "inputs": [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for path in (BATCH, PENCIL, JACOBIAN, MINIMAL)
    ],
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "explicit generic birational maps in both directions over GF(19^2)(V)",
            "literal generic function-field replay, not only invariant agreement",
            "explicit composition with the minimal short Weierstrass gauge",
        ],
        "not_proved": [
            "transported old-component and zero-section marking",
            "Frobenius-invariant normalization or a second-prime alignment",
            "characteristic-zero Mordell--Weil rank",
        ],
    },
    "reproduce": "sage -python elkies-k3/scripts/interpolate_q80_third_q12_maps_mod19_quadratic.sage",
}
OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
forward_degrees = [
    record["degrees_V_numerator_denominator"]
    for records in forward_records.values()
    for record in records
]
print(
    f"Q80THIRDQ12MAPS|training={len(training)}|heldout={len(held_out)}|"
    f"forward_V_degrees={forward_degrees}|"
    f"inverse_bounds={[record['weighted_bound'] for record in inverse_records.values()]}|"
    "generic_forward=PASS|generic_inverse=PASS|"
    "status=PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_MOD19_QUADRATIC",
    flush=True,
)
