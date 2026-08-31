#!/usr/bin/env sage
"""Construct the q52 horizontal/root seeds by fibrewise Abel words mod 131.

Each of the selected old 3A1 sections is a degree-15--25 multisection of the
physical q114 4A1 pencil and has an enormous individual Abel image.  Reduce
the multisections fibrewise and add the integral q52 words before rational
interpolation.  This preserves the cancellations giving P.O=15,1,1,3 and
avoids constructing any huge intermediate section.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SOURCE_SURFACE = LOCAL / "fixed-reverse-3a1-rr-qq.json"
SOURCE_CURVES = LOCAL / "fixed-reverse-4a1-horizontal-from-3a1-qq.json"
CURRENT_RR = LOCAL / "fixed-reverse-4a1-rr-qq.json"
CURRENT_POINTING = LOCAL / "fixed-reverse-4a1-pointing-qq.json"
CURRENT_AUDIT = LOCAL / "fixed-reverse-4a1-physical-nef-audit.json"
Q52_AUDIT = LOCAL / "fixed-reverse-5a1-physical-nef-audit.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=131)
parser.add_argument("--tau", type=int, default=7)
parser.add_argument("--interpolate", action="store_true")
parser.add_argument(
    "--good-sample-target", type=int, default=110,
    help="stop modular interpolation after this many good fibres",
)
parser.add_argument(
    "--sample-attempt-limit", type=int, default=10000,
    help="maximum consecutive integer specializations to try",
)
parser.add_argument("--output", type=Path)
args = parser.parse_args()


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


started = time.monotonic()
source_surface = read_json(SOURCE_SURFACE)
source_curves = read_json(SOURCE_CURVES)
current_rr = read_json(CURRENT_RR)
current_pointing = read_json(CURRENT_POINTING)
current_audit = read_json(CURRENT_AUDIT)
q52_audit = read_json(Q52_AUDIT)
manifest = read_json(MANIFEST)
assert current_rr["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_RR_JACOBIAN"
assert current_pointing["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_POINTING"
assert q52_audit["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_5A1_PHYSICAL_NEF"

p = ZZ(args.prime)
F = GF(p)
R = PolynomialRing(F, "s")
s = R.gen()
K = R.fraction_field()
T = PolynomialRing(F, "t")
t = T.gen()
KT = T.fraction_field()
U = PolynomialRing(T, "s")
ss = U.gen()


def red(value):
    value = QQ(value)
    return F(value.numerator()) / F(value.denominator())


def poly(values):
    return R([red(value) for value in values])


def tpoly(values):
    return T([red(value) for value in values])


def rf_record(record, ring=R, field=K):
    return field(ring([red(value) for value in record["numerator_coefficients_low_to_high"]])) / field(
        ring([red(value) for value in record["denominator_coefficients_low_to_high"]])
    )


def point_record(record):
    x = K(poly(record["x_numerator_coefficients_low_to_high"])) / K(
        poly(record["x_denominator_coefficients_low_to_high"])
    )
    y = K(poly(record["y_numerator_coefficients_low_to_high"])) / K(
        poly(record["y_denominator_coefficients_low_to_high"])
    )
    return x, y


A3 = poly(source_surface["child"]["minimal_A_coefficients_low_to_high"])
B3 = poly(source_surface["child"]["minimal_B_coefficients_low_to_high"])
horizontal = source_curves["section"]
X = poly(horizontal["x_numerator_coefficients_low_to_high"])
Y = poly(horizontal["y_numerator_coefficients_low_to_high"])
Z = poly(horizontal["Z_coefficients_low_to_high"])
Hx = K(X) / K(Z ** 2)
Hy = K(Y) / K(Z ** 3)
basis = current_rr["smooth_RR"]["basis_pairs"]
AA0 = poly(basis[0]["AA_coefficients_low_to_high"])
BB0 = poly(basis[0]["BB_coefficients_low_to_high"])
AA1 = poly(basis[1]["AA_coefficients_low_to_high"])
BB1 = poly(basis[1]["BB_coefficients_low_to_high"])
quartic = sum((
    tpoly(values) * ss ** degree
    for degree, values in enumerate(
        current_rr["binary_quartic"]["coefficients_in_old_u_low_to_high"]
    )
), U.zero())
square_factor = sum((
    tpoly(values) * ss ** degree
    for degree, values in enumerate(
        current_rr["binary_quartic"]["square_factor_coefficients_in_old_u_low_to_high"]
    )
), U.zero())
alpha = red(QQ(current_pointing["fixed_zero"]["old_I2_support"]))
q_origin = rf_record(current_pointing["fixed_zero"]["quartic_ordinate"], T, KT)
A4 = tpoly(current_rr["child"]["minimal_A_coefficients_low_to_high"])
B4 = tpoly(current_rr["child"]["minimal_B_coefficients_low_to_high"])
I2_factor = next(
    tpoly(record["factor_coefficients_low_to_high"])
    for record in current_rr["child"]["finite_fibres"]
    if record["kodaira"] == "I2"
)
I2_supports = I2_factor.roots(multiplicities=False)
assert len(I2_supports) == 4


def eval_t_at_base(value, base):
    value = KT(value)

    def evaluate(poly_value):
        answer = K.zero()
        for coefficient in reversed(T(poly_value).list()):
            answer = answer * base + K(coefficient)
        return answer

    return evaluate(value.numerator()) / evaluate(value.denominator())


def eval_bivariate(poly_value, base):
    answer = K.zero()
    for coefficient in reversed(U(poly_value).list()):
        answer = answer * K(s) + eval_t_at_base(coefficient, base)
    return answer


e = KT(U(quartic)(U(alpha)))
d = KT(U(quartic.derivative())(U(alpha)))
c = KT(U(quartic.derivative(2))(U(alpha))) / 2
b = KT(U(quartic.derivative(3))(U(alpha))) / 6
assert e == q_origin ** 2
a1 = d / q_origin
a2 = c - d ** 2 / (4 * q_origin ** 2)
a3 = 2 * q_origin * b
b2 = a1 ** 2 + 4 * a2


def parametrized_child_curve(P):
    px, py = P
    assert py ** 2 == px ** 3 + K(A3) * px + K(B3)
    slope = (py + Hy) / (px - Hx)
    base = -(K(AA0) + K(BB0 * Z) * slope) / (
        K(AA1) + K(BB1 * Z) * slope
    )
    bb_value = K(BB0) + base * K(BB1)
    ordinate = (
        bb_value ** 2 * (2 * px + Hx - slope ** 2)
        / eval_bivariate(square_factor, base)
    )
    assert ordinate ** 2 == eval_bivariate(quartic, base)
    q_value = eval_t_at_base(q_origin, base)
    d_value = eval_t_at_base(d, base)
    c_value = eval_t_at_base(c, base)
    a1_value = eval_t_at_base(a1, base)
    a3_value = eval_t_at_base(a3, base)
    b2_value = eval_t_at_base(b2, base)
    relative = K(s - alpha)
    x_general = (
        2 * q_value * (ordinate + q_value) + d_value * relative
    ) / relative ** 2
    y_general = (
        4 * q_value ** 2 * (ordinate + q_value)
        + 2 * q_value * d_value * relative
        + (2 * q_value * c_value - d_value ** 2 / (2 * q_value)) * relative ** 2
    ) / relative ** 3
    x_child = K(9 * (x_general + b2_value / 12))
    y_child = K(27 * (y_general + (a1_value * x_general + a3_value) / 2))
    assert y_child ** 2 == x_child ** 3 + eval_t_at_base(A4, base) * x_child + eval_t_at_base(B4, base)
    return base, x_child, y_child


# Lattice words in the fourteen stored source curves.  These are fixed
# geometric curves, so express them in the reflected physical 4A1 basis;
# merely applying the raw q114 transition would use the wrong marking.
physical_basis = matrix(
    ZZ,
    current_audit["full_marked_transport"][
        "physical_4A1_basis_in_3A1_coordinates"
    ],
)
classes = [
    vector(ZZ, record["class_in_3A1_coordinates"]) * physical_basis.inverse()
    for record in source_curves["candidate_construction"]
]
tail_generators = matrix(ZZ, [list(value[6:]) for value in classes])
smith, left_transform, right_transform = tail_generators.smith_form(
    transformation=True
)
assert left_transform * tail_generators * right_transform == smith
assert smith[:13, :] == matrix.identity(ZZ, 13)
kernel_word = vector(ZZ, left_transform.row(13))
assert kernel_word * tail_generators == 0
target_section = vector(ZZ, [
    16, 1, -3, -7, 4, 0, 2, -3, -2, 0, 0, 0, 1, -1, 3, 0, 1, 0, 0,
])
root_classes = []
for values in q52_audit["full_marked_transport"]["physical_5A1_root_classes_in_4A1_coordinates"]:
    root = vector(ZZ, values)
    if abs(root[1]) == 1:
        root_classes.append(root if root[1] == 1 else -root)
targets = [target_section] + root_classes
words = []
for target in targets:
    transformed_target = vector(ZZ, target[6:]) * right_transform
    word = vector(ZZ, list(transformed_target) + [0]) * left_transform
    word = min(
        (word + multiple * kernel_word for multiple in range(-100, 101)),
        key=lambda value: sum(abs(entry) for entry in value),
    )
    assert word * tail_generators == vector(ZZ, target[6:])
    words.append(word)
needed_positions = sorted({
    position for word in words for position, coefficient in enumerate(word) if coefficient
})
needed_source_indices = needed_positions
source_points = [
    point_record(source_curves["exact_3A1_MW_basis"][index]["section"])
    for index in needed_source_indices
]
parametric_curves = [parametrized_child_curve(point) for point in source_points]


def reduce_mod_H(value, H):
    value = K(value)
    numerator = R(value.numerator())
    denominator = R(value.denominator())
    if denominator.gcd(H).degree() != 0:
        raise ZeroDivisionError("multisection denominator is not invertible")
    return (numerator * denominator.inverse_mod(H)) % H


def newton_power_sums(polynomial):
    degree = polynomial.degree()
    assert polynomial[degree] == 1
    sums = [F(degree)]
    for order in range(1, degree):
        total = F(order) * polynomial[degree - order]
        for index in range(1, order):
            total += polynomial[degree - index] * sums[order - index]
        sums.append(-total)
    return sums


def abel_trace_at(curve, tau, E_tau):
    base, x_curve, y_curve = curve
    H = R(base.numerator() - tau * base.denominator())
    degree = H.degree()
    generic_degree = max(
        R(base.numerator()).degree(), R(base.denominator()).degree()
    )
    if degree != generic_degree:
        raise ArithmeticError("multisection has a point at source infinity")
    H = H.monic()
    if H.gcd(H.derivative()).degree() != 0:
        raise ArithmeticError("multisection fibre is not etale")
    xA = reduce_mod_H(x_curve, H)
    yA = reduce_mod_H(y_curve, H)
    Atau = F(A4(tau))
    Btau = F(B4(tau))
    if (yA ** 2 - xA ** 3 - Atau * xA - Btau) % H:
        raise ArithmeticError("parametric curve misses the child fibre")
    pole_bound = degree + 1
    x_count = pole_bound // 2 + 1
    y_count = (pole_bound - 3) // 2 + 1 if pole_bound >= 3 else 0
    xp = [R.one()]
    for unused in range(max(x_count, y_count) - 1):
        xp.append((xp[-1] * xA) % H)
    columns = xp[:x_count] + [(yA * xp[index]) % H for index in range(y_count)]
    assert len(columns) == degree + 1
    evaluation = matrix(F, degree, degree + 1, lambda row, column: columns[column][row])
    kernel = evaluation.right_kernel().basis_matrix()
    if kernel.nrows() != 1:
        raise ArithmeticError("Abel kernel is not one-dimensional")
    relation = kernel[0]
    XR = PolynomialRing(F, "X")
    XX = XR.gen()
    Afun = sum(relation[index] * XX ** index for index in range(x_count))
    Bfun = sum(relation[x_count + index] * XX ** index for index in range(y_count))
    intersection = Afun ** 2 - (XX ** 3 + Atau * XX + Btau) * Bfun ** 2
    if intersection.degree() != degree + 1:
        raise ArithmeticError("Abel residual intersection has the wrong degree")
    root_sum = -intersection[degree] / intersection[degree + 1]
    power_sums = newton_power_sums(H)
    divisor_x_sum = sum(xA[index] * power_sums[index] for index in range(degree))
    residual_x = root_sum - divisor_x_sum
    if not Bfun(residual_x):
        raise ArithmeticError("Abel residual has vanishing B coordinate")
    residual_y = -Afun(residual_x) / Bfun(residual_x)
    return -E_tau(residual_x, residual_y), degree


def multiply(P, coefficient):
    coefficient = ZZ(coefficient)
    return coefficient * P


def word_points_at(tau):
    tau = F(tau)
    Atau, Btau = F(A4(tau)), F(B4(tau))
    E_tau = EllipticCurve(F, [0, 0, 0, Atau, Btau])
    if not E_tau.discriminant():
        raise ArithmeticError("singular child fibre")
    traces = []
    degrees = []
    for curve in parametric_curves:
        trace, degree = abel_trace_at(curve, tau, E_tau)
        traces.append(trace)
        degrees.append(degree)
    kernel_check = E_tau.zero()
    for trace_position, source_position in enumerate(needed_positions):
        kernel_check += multiply(traces[trace_position], kernel_word[source_position])
    if not kernel_check.is_zero():
        raise ArithmeticError("physical MW lattice relation failed on the fibre")
    results = []
    for word in words:
        answer = E_tau.zero()
        for trace_position, source_position in enumerate(needed_positions):
            answer += multiply(traces[trace_position], word[source_position])
        if answer.is_zero():
            raise ArithmeticError("q52 word specialized to zero")
        results.append(tuple(map(F, answer.xy())))
    return results, degrees


def interpolate_rational(samples, target_index, coordinate_index, numerator_degree, denominator_degree):
    needed = numerator_degree + denominator_degree + 1
    if len(samples) <= needed:
        raise ArithmeticError("not enough samples for interpolation and holdout")
    training = samples[:needed]
    rows = []
    for tau, points in training:
        value = points[target_index][coordinate_index]
        rows.append(
            [tau ** degree for degree in range(numerator_degree + 1)]
            + [-value * tau ** degree for degree in range(denominator_degree + 1)]
        )
    kernel = matrix(F, rows).right_kernel().basis_matrix()
    if kernel.nrows() != 1:
        raise ArithmeticError("interpolation kernel is not one-dimensional")
    relation = kernel[0]
    numerator = T(list(relation[:numerator_degree + 1]))
    denominator = T(list(relation[numerator_degree + 1:]))
    if denominator.degree() != denominator_degree:
        raise ArithmeticError("interpolation denominator dropped degree")
    for tau, points in samples[needed:]:
        if denominator(tau) and numerator(tau) / denominator(tau) != points[target_index][coordinate_index]:
            raise ArithmeticError(
                "interpolation holdout failed for target {} coordinate {} at {}"
                .format(target_index, coordinate_index, tau)
            )
    return numerator, denominator, needed, len(samples) - needed


def homogeneous_seed_and_rank(x, y, pole):
    x_denominator = T(x.denominator())
    monic_denominator = x_denominator / x_denominator.leading_coefficient()
    if not monic_denominator.is_square():
        raise ArithmeticError("interpolated x denominator is not a square")
    Z = T(monic_denominator.sqrt())
    X = T(x * KT(Z ** 2))
    Y = T(y * KT(Z ** 3))
    assert Z.is_monic() and Z.degree() == pole
    assert Y ** 2 == X ** 3 + A4 * X * Z ** 4 + B4 * Z ** 6
    derivatives = []
    derivatives.extend(
        -(3 * X ** 2 + A4 * Z ** 4) * t ** degree
        for degree in range(2 * pole + 5)
    )
    derivatives.extend(
        2 * Y * t ** degree for degree in range(3 * pole + 7)
    )
    derivatives.extend(
        -(4 * A4 * X * Z ** 3 + 6 * B4 * Z ** 5) * t ** degree
        for degree in range(pole)
    )
    equation_degree = 6 * pole + 12
    rows = [
        [
            derivative[degree] if degree <= derivative.degree() else F.zero()
            for derivative in derivatives
        ]
        for degree in range(equation_degree + 1)
    ]
    coefficient_jacobian = matrix(F, rows)
    coefficient_rank = int(coefficient_jacobian.rank())
    node_hits = []
    for support in I2_supports:
        node = -3 * B4(support) / (2 * A4(support))
        if Y(support) or X(support) != node * Z(support) ** 2:
            continue
        node_hits.append((support, node))
        rows.append(
            [support ** degree for degree in range(2 * pole + 5)]
            + [F.zero()] * (3 * pole + 7)
            + [
                -2 * node * Z(support) * support ** degree
                for degree in range(pole)
            ]
        )
        rows.append(
            [F.zero()] * (2 * pole + 5)
            + [support ** degree for degree in range(3 * pole + 7)]
            + [F.zero()] * pole
        )
    jacobian = matrix(F, rows)
    rank = int(jacobian.rank())
    pivot_rows = list(map(int, jacobian.transpose().pivots()))
    determinant = None
    if rank == len(derivatives):
        determinant = int(
            matrix(F, [jacobian.row(row) for row in pivot_rows]).det()
        )
        assert determinant
    return (
        X, Y, Z, coefficient_rank, rank, pivot_rows, determinant,
        [(int(support), int(node)) for support, node in node_hits],
    )


expected_poles = [15, 1, 1, 3]
if args.interpolate:
    samples = []
    degree_fingerprint = None
    attempt_count = min(int(p), int(args.sample_attempt_limit))
    for tau_integer in range(attempt_count):
        tau_value = F(tau_integer)
        try:
            points, degrees = word_points_at(tau_value)
            samples.append((F(tau_value), points))
            degree_fingerprint = degrees
            if len(samples) >= int(args.good_sample_target):
                break
        except (ArithmeticError, ZeroDivisionError, ValueError):
            continue
    section_records = []
    for target_index, pole in enumerate(expected_poles):
        xn, xd, xtrain, xholdout = interpolate_rational(
            samples, target_index, 0, 2 * pole + 4, 2 * pole
        )
        yn, yd, ytrain, yholdout = interpolate_rational(
            samples, target_index, 1, 3 * pole + 6, 3 * pole
        )
        x = KT(xn) / KT(xd)
        y = KT(yn) / KT(yd)
        assert y ** 2 == x ** 3 + KT(A4) * x + KT(B4)
        (
            X, Y, Z, coefficient_rank, rank, pivot_rows, determinant,
            node_hits,
        ) = homogeneous_seed_and_rank(x, y, pole)
        section_records.append({
            "P_dot_O": pole,
            "x_numerator_coefficients_low_to_high": list(map(int, xn.list())),
            "x_denominator_coefficients_low_to_high": list(map(int, xd.list())),
            "y_numerator_coefficients_low_to_high": list(map(int, yn.list())),
            "y_denominator_coefficients_low_to_high": list(map(int, yd.list())),
            "degrees_x_num_den_y_num_den": [
                int(xn.degree()), int(xd.degree()), int(yn.degree()), int(yd.degree())
            ],
            "training_holdout": [xtrain, xholdout, ytrain, yholdout],
            "X_coefficients_low_to_high": list(map(int, X.list())),
            "Y_coefficients_low_to_high": list(map(int, Y.list())),
            "Z_coefficients_low_to_high": list(map(int, Z.list())),
            "coefficient_jacobian_rank_before_node_resolution": coefficient_rank,
            "resolved_jacobian_rank": rank,
            "coefficient_variable_count": 6 * pole + 12,
            "I2_node_hits_support_and_node": node_hits,
            "selected_independent_equation_rows": pivot_rows,
            "selected_jacobian_determinant": determinant,
        })
    status = (
        "PASS_MODP_FIXED_REVERSE_5A1_ABEL_WORD_REGULAR_SEEDS"
        if all(
            record["resolved_jacobian_rank"]
            == record["coefficient_variable_count"]
            for record in section_records
        )
        else "PASS_MODP_FIXED_REVERSE_5A1_ABEL_WORD_SINGULAR_SEEDS"
    )
    body = {
        "sections": section_records,
        "good_fibre_count": len(samples),
        "sample_attempt_limit": attempt_count,
        "good_sample_target": int(args.good_sample_target),
        "multisection_degrees": list(map(int, degree_fingerprint)),
    }
else:
    points, degrees = word_points_at(F(args.tau))
    status = "PASS_MODP_FIXED_REVERSE_5A1_ABEL_WORD_SMOKE"
    body = {
        "tau": int(F(args.tau)),
        "points": [[int(x), int(y)] for x, y in points],
        "multisection_degrees": list(map(int, degrees)),
    }

payload = {
    "schema": "elkies-k3.fixed-reverse-5a1-abel-word-modp.v1",
    "status": status,
    "prime": int(p),
    "source_generator_count": len(classes),
    "needed_source_basis_indices": list(map(int, needed_source_indices)),
    "words": [list(map(int, word)) for word in words],
    **body,
    "method": {
        "fibrewise_abel_reduction_before_word_sum": True,
        "individual_huge_abel_sections_constructed": False,
        "groebner_or_surface_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "Exact finite-field construction and equation verification. Characteristic-zero "
        "lifting, exact QQ section identities, RR compilation and pointing remain separate."
    ),
    "inputs": {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (
            SOURCE_SURFACE, SOURCE_CURVES, CURRENT_RR, CURRENT_POINTING,
            Q52_AUDIT, MANIFEST,
        )
    },
}
default_name = (
    "fixed-reverse-5a1-abel-word-seeds-mod{}.json".format(p)
    if args.interpolate else
    "fixed-reverse-5a1-abel-word-smoke-mod{}.json".format(p)
)
output = args.output or (LOCAL / default_name)
output = output if output.is_absolute() else ROOT / output
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE5A1ABEL|p={}|interpolate={}|needed={}|degrees={}|seconds={:.3f}|status={}|output={}".format(
        p, int(args.interpolate), len(needed_source_indices), body["multisection_degrees"],
        payload["method"]["runtime_seconds"], status, output,
    ),
    flush=True,
)
