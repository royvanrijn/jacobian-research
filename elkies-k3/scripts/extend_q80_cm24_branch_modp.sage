#!/usr/bin/env sage
"""Continue a characteristic-zero q=80 CM24 tangent after split reduction.

The exact CM24 marked system is defined over Q(sqrt(-3),sqrt(-6)).  At a
prime where both radicals split, this script reduces that *same* system,
selects a tangent by its rational surface coordinates (D',P'), and performs
the normalized formal lift P=P0+h.  Its main use is to distinguish genuine
reductions of the characteristic-zero branch from phenomena special to 7.

Exported relation spaces are canonical reduced-row-echelon bases on the
declared centered monomial ordering.  They are finite modular evidence, not
characteristic-zero identities.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing, PowerSeriesRing, QQ, ZZ, binomial, vector


parser = argparse.ArgumentParser()
parser.add_argument(
    "--prime", type=int, default=31,
    help="split prime, or 0 for the exact CM number field",
)
parser.add_argument("--order", type=int, default=214)
parser.add_argument("--slope", choices=("8/87", "1/12"), default="1/12")
parser.add_argument("--relation-degree", type=int, default=6)
parser.add_argument("--relation-basis-output")
parser.add_argument("--jet-output")
parser.add_argument("--parameter-input")
parser.add_argument("--algebraic-coordinate")
parser.add_argument("--algebraic-coordinate-degree", type=int, default=2)
parser.add_argument("--algebraic-max-parameter-degree", type=int, default=40)
parser.add_argument("--algebraic-validation", type=int, default=8)
arguments = parser.parse_args()
if arguments.order < 4:
    parser.error("--order must be at least 4")
if arguments.relation_degree < 0:
    parser.error("--relation-degree must be nonnegative")

# This is the canonical exact boundary/system certificate.  Everything below
# is obtained by reducing its constants, not by rebuilding a separate mod-p
# ansatz.
load("elkies-k3/scripts/verify_q80_cm24_rational_model.sage")

prime = ZZ(arguments.prime)
exact_mode = prime == 0
if exact_mode:
    field = composite
    r3 = s3
    r6 = s6
else:
    if not prime.is_prime() or prime in (2, 3, 13):
        raise ValueError("choose a good prime other than 2, 3, or 13")
    field = GF(prime)
    roots_three = field(-3).sqrt(all=True)
    roots_six = field(-6).sqrt(all=True)
    if len(roots_three) != 2 or len(roots_six) != 2:
        raise ValueError("both -3 and -6 must split at the selected prime")
    r3 = min(roots_three, key=ZZ)
    r6 = min(roots_six, key=ZZ)


def reduce_value(value):
    if exact_mode:
        return composite(value)
    coordinates = biquadratic_coordinates(composite(value))
    return sum(
        (field(coordinates[0]), field(coordinates[1])*r3,
         field(coordinates[2])*r6, field(coordinates[3])*r3*r6),
        field.zero(),
    )


# Reuse the exact polynomial evaluator with a finite base ring.  Its explicit
# rational constants coerce into this field; all CM constants enter through
# the reduced seed.
Jet2.base_ring = field
seed = vector(field, [reduce_value(value) for value in cone_seed_values])
Jet2.precision = 2
assert not any(
    value.constant for value in evaluate_cone_system(list(map(Jet2, seed)))
)

jacobian_columns = []
for column in range(len(seed)):
    inputs = [
        Jet2(value, 1 if index == column else 0)
        for index, value in enumerate(seed)
    ]
    jacobian_columns.append(
        vector(field, [value.linear for value in evaluate_cone_system(inputs)])
    )
jacobian = Matrix(field, jacobian_columns).transpose()
resolved = jacobian.matrix_from_columns(cone_active_columns)
kernel = resolved.right_kernel_matrix()
left_kernel = resolved.left_kernel_matrix()
if resolved.rank() != 37 or kernel.nrows() != 2:
    raise ValueError(
        f"bad reduction: rank={resolved.rank()} nullity={kernel.nrows()}"
    )

active_seed = vector(field, [seed[column] for column in cone_active_columns])
active_names = tuple(cone_active_names)
d_column = active_names.index("D")
p_column = active_names.index("P")
slope = field(QQ(arguments.slope))

# Find the unique tangent vector with D'=slope and P'=1.  This avoids any
# dependence on the basis chosen for the two-dimensional tangent kernel.
surface_projection = Matrix(
    field,
    [[kernel[row, column] for row in range(2)] for column in (d_column, p_column)],
)
if surface_projection.rank() != 2:
    raise ValueError("surface tangent projection degenerated")
tangent_coordinates = surface_projection.solve_right(vector(field, [slope, 1]))
first_order = tangent_coordinates * kernel
assert first_order[d_column] == slope and first_order[p_column] == 1


series_ring = PowerSeriesRing(field, "h", default_prec=arguments.order+2)


def series_poly_add(left, right):
    length = max(len(left), len(right))
    return [
        (left[index] if index < len(left) else series_ring.zero())
        +(right[index] if index < len(right) else series_ring.zero())
        for index in range(length)
    ]


def series_poly_neg(value):
    return [-coefficient for coefficient in value]


def series_poly_mul(left, right, truncate=None):
    length = len(left)+len(right)-1
    if truncate is not None:
        length = min(length, truncate)
    result = [series_ring.zero() for _ in range(length)]
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            if left_index+right_index < length:
                result[left_index+right_index] += left_value*right_value
    return result


def series_poly_pow(value, exponent, truncate=None):
    result = [series_ring.one()]
    for _ in range(exponent):
        result = series_poly_mul(result, value, truncate=truncate)
    return result


def series_poly_scale(value, scalar):
    scalar = series_ring(scalar)
    return [scalar*coefficient for coefficient in value]


def series_poly_subtract(first, *others):
    result = first
    for other in others:
        result = series_poly_add(result, series_poly_neg(other))
    return result


def evaluate_series_system(values):
    """Power-series version of the pinned exact cone evaluator."""
    D, P, Q, E = values[:4]
    remaining = values[4:]
    X1 = list(remaining[0:5])
    Y1 = list(remaining[5:12])
    X2 = list(remaining[12:17])
    Y2 = list(remaining[17:24])
    LAM = remaining[24]
    N = list(remaining[25:32])
    M = list(remaining[32:42])
    R = -3*D**2+3-P-Q
    A_series = [series_ring.zero(), series_ring.zero(), series_ring(-3), P, Q, R]
    A_at_one = [
        sum(
            (series_ring(binomial(index, order))*A_series[index]
             for index in range(order, 6)),
            series_ring.zero(),
        )
        for order in range(4)
    ]
    denominator = -3*D**2
    u_series = list(A_at_one)
    u_series[0] += 3*D**2
    u_series = [coefficient/denominator for coefficient in u_series]
    branch_series = [series_ring.one(), series_ring.zero(), series_ring.zero(), series_ring.zero()]
    branch_series = series_poly_add(
        branch_series, series_poly_scale(u_series, field(3)/2)
    )
    branch_series = series_poly_add(
        branch_series,
        series_poly_scale(series_poly_pow(u_series, 2, truncate=4), field(3)/8),
    )
    branch_series = series_poly_add(
        branch_series,
        series_poly_scale(series_poly_pow(u_series, 3, truncate=4), -field(1)/16),
    )
    branch_series = series_poly_scale(branch_series, 2*D**3)
    fixed_series = [
        series_ring(2*binomial(3, order))+E*binomial(8, order)
        for order in range(4)
    ]
    difference = [branch_series[index]-fixed_series[index] for index in range(4)]
    middle_B = [
        sum(
            (series_ring(field(jet_matrix_inverse[row, column]))*difference[column]
             for column in range(4)),
            series_ring.zero(),
        )
        for row in range(4)
    ]
    B_series = [series_ring.zero()]*3+[series_ring(2)]+middle_B+[E]
    H = [-LAM, series_ring.one()]
    residuals = (
        series_poly_subtract(
            series_poly_pow(Y1, 2), series_poly_pow(X1, 3),
            series_poly_mul(A_series, X1), B_series,
        ),
        series_poly_subtract(
            series_poly_pow(Y2, 2), series_poly_pow(X2, 3),
            series_poly_mul(A_series, X2), B_series,
        ),
        series_poly_subtract(
            series_poly_pow(M, 2), series_poly_pow(N, 3),
            series_poly_mul(series_poly_mul(A_series, N), series_poly_pow(H, 4)),
            series_poly_mul(B_series, series_poly_pow(H, 6)),
        ),
    )
    equations = []
    for residual, expected_degree in zip(residuals, (12, 12, 18)):
        residual += [series_ring.zero()]*(expected_degree+1-len(residual))
        equations.extend(residual[:expected_degree+1])
    value_at_one = lambda value: sum(value, series_ring.zero())
    derivative_at_one = lambda value: sum(
        (series_ring(index)*value[index] for index in range(1, len(value))),
        series_ring.zero(),
    )
    X1_one = value_at_one(X1)
    Y1_one = value_at_one(Y1)
    N_one = value_at_one(N)
    M_one = value_at_one(M)
    A_prime = derivative_at_one(A_series)
    X1_prime = derivative_at_one(X1)
    Y1_prime = derivative_at_one(Y1)
    H_one = 1-LAM
    M_affine_prime = derivative_at_one(M)/H_one**3-3*M_one/H_one**4
    N_affine_prime = derivative_at_one(N)/H_one**2-2*N_one/H_one**3
    equations.extend((
        X1_one-D,
        Y1_one,
        N_one-D*H_one**2,
        M_one,
        X1[1]-1,
        (6*D*Y1_prime)**2-3*D*(6*D*X1_prime+A_prime)**2,
        (6*D*M_affine_prime)**2-3*D*(6*D*N_affine_prime+A_prime)**2,
    ))
    assert len(equations) == 52
    return equations


def residual_values(coefficients, precision):
    active_lookup = {
        column: index for index, column in enumerate(cone_active_columns)
    }
    inputs = []
    for column, seed_value in enumerate(seed):
        if column in active_lookup:
            active_index = active_lookup[column]
            series = [row[active_index] for row in coefficients]
        else:
            series = [seed_value]
        inputs.append(series_ring(series).add_bigoh(precision))
    return evaluate_series_system(inputs)


def residual_coefficient(coefficients, order):
    return vector(
        field,
        [value[order] for value in residual_values(coefficients, order+1)],
    )


coefficients = [active_seed, vector(field, first_order)]
for order in range(2, arguments.order):
    residual = residual_coefficient(coefficients, order)
    assert not any(left_kernel*residual)
    particular = resolved.solve_right(-residual)

    # At higher order P is fixed to zero.  This cuts the tangent kernel to a
    # line.  The look-ahead obstruction at order+1 fixes that last scalar.
    normalization_row = Matrix(
        field, [[kernel[row, p_column] for row in range(2)]]
    )
    normalization_kernel = normalization_row.right_kernel_matrix()
    assert normalization_kernel.nrows() == 1
    free_tangent = normalization_kernel.row(0)*kernel
    assert not free_tangent[p_column]
    correction = particular
    if correction[p_column]:
        coordinate_fix = normalization_row.solve_right(
            vector(field, [-correction[p_column]])
        )

        correction += coordinate_fix*kernel
    assert not correction[p_column]

    trial0 = coefficients+[vector(field, correction)]
    trial1 = coefficients+[vector(field, correction+free_tangent)]
    obstruction0 = left_kernel*residual_coefficient(trial0, order+1)
    obstruction_step = (
        left_kernel*residual_coefficient(trial1, order+1)-obstruction0
    )
    forced = {
        -constant/linear
        for constant, linear in zip(obstruction0, obstruction_step)
        if linear
    }
    assert all(
        linear or not constant
        for constant, linear in zip(obstruction0, obstruction_step)
    )
    assert len(forced) == 1, (order, forced)
    scalar = forced.pop()
    coefficients.append(vector(field, correction+scalar*free_tangent))
    assert not any(
        left_kernel*residual_coefficient(coefficients, order+1)
    )
    if (order+1) % 10 == 0 or order+1 == arguments.order:
        print(
            f"Q80CM24MODP|prime={'QQCM' if exact_mode else prime}|slope={arguments.slope}|"
            f"stage=lift|order={order+1}|unique=1",
            flush=True,
        )

if arguments.algebraic_coordinate:
    if exact_mode or not arguments.parameter_input:
        raise ValueError("algebraic coordinate fitting needs a finite prime and --parameter-input")
    coordinate_name = arguments.algebraic_coordinate
    if coordinate_name not in active_names:
        raise ValueError("unknown active algebraic coordinate")
    parameter_artifact = json.loads(Path(arguments.parameter_input).read_text())
    if parameter_artifact.get("schema") != "q80-cm24-formal-branch-parameter-v1":
        raise ValueError("unexpected surface parameter artifact")
    expected_slope = "5:1" if arguments.slope == "8/87" and prime == 7 else None
    if expected_slope and parameter_artifact.get("slope_mod_7") != expected_slope:
        raise ValueError("surface parameter belongs to another branch")
    parameter_series_ring = PowerSeriesRing(field, "t", default_prec=arguments.order)
    t_series = parameter_series_ring.gen()
    parameter_polynomials = PolynomialRing(field, "t")
    parameter_functions = parameter_polynomials.fraction_field()
    centered_p = parameter_functions(parameter_artifact["functions"]["P"]["value"])
    h_of_t = (
        parameter_series_ring(centered_p.numerator()(t_series))
        /parameter_series_ring(centered_p.denominator()(t_series))
    ).add_bigoh(arguments.order)
    if h_of_t.valuation() != 1:
        raise ValueError("surface parameter is not local at CM24")
    h_coefficients = tuple(
        row[active_names.index(coordinate_name)] for row in coefficients
    )
    coordinate_of_t = parameter_series_ring.zero()
    power = parameter_series_ring.one()
    for coefficient in h_coefficients:
        coordinate_of_t += coefficient*power
        power = (power*h_of_t).add_bigoh(arguments.order)
    coordinate_of_t = coordinate_of_t.add_bigoh(arguments.order)
    coordinate_degree = arguments.algebraic_coordinate_degree
    validation = arguments.algebraic_validation
    fitting_order = arguments.order-validation
    relation_ring = PolynomialRing(field, names=("t", coordinate_name))
    t_variable, coordinate_variable = relation_ring.gens()
    coordinate_powers = [parameter_series_ring.one()]
    for _ in range(coordinate_degree):
        coordinate_powers.append(
            (coordinate_powers[-1]*coordinate_of_t).add_bigoh(arguments.order)
        )
    relation = None
    selected_parameter_degree = None
    for parameter_degree in range(arguments.algebraic_max_parameter_degree+1):
        monomials = tuple(
            t_variable**parameter_exponent*coordinate_variable**coordinate_exponent
            for coordinate_exponent in range(coordinate_degree+1)
            for parameter_exponent in range(parameter_degree+1)
        )
        if len(monomials) >= fitting_order:
            continue
        columns = []
        for monomial in monomials:
            parameter_exponent, coordinate_exponent = monomial.exponents()[0]
            value = (
                t_series**parameter_exponent*coordinate_powers[coordinate_exponent]
            ).add_bigoh(arguments.order)
            columns.append(vector(field, [value[index] for index in range(fitting_order)]))
        kernel_relation = Matrix(field, columns).transpose().right_kernel_matrix()
        if kernel_relation.nrows() != 1:
            continue
        candidate = sum(
            coefficient*monomial
            for coefficient, monomial in zip(kernel_relation.row(0), monomials)
        )
        if candidate.degree(coordinate_variable) != coordinate_degree:
            continue
        residual = parameter_series_ring.zero()
        for exponents, coefficient in candidate.dict().items():
            residual += (
                coefficient*t_series**exponents[0]*coordinate_powers[exponents[1]]
            )
        if residual.add_bigoh(arguments.order):
            continue
        relation = candidate
        selected_parameter_degree = parameter_degree
        break
    print(
        f"Q80CM24MODP|prime={prime}|slope={arguments.slope}|"
        f"stage=algebraic_coordinate|coordinate={coordinate_name}|"
        f"coordinate_degree={coordinate_degree}|"
        f"parameter_degree={selected_parameter_degree}|validation={validation}|"
        f"relation={relation}|status={'PASS' if relation is not None else 'NO_RELATION'}",
        flush=True,
    )

surface_names = ("D", "P", "Q", "E")
surface_columns = tuple(active_names.index(name) for name in surface_names)
centers = tuple(field(value) for value in (QQ(-1)/2, QQ(9)/4, QQ(-9)/4, QQ(-27)/32))
surface_series = []
for column, center in zip(surface_columns, centers):
    values = [row[column] for row in coefficients]
    values[0] -= center
    if exact_mode:
        # The unmarked surface branch is rational although its selected
        # sections live in the CM compositum.
        values = [QQ(value) for value in values]
    surface_series.append(tuple(values))


def truncated_product(left, right, length):
    base = QQ if exact_mode else field
    answer = [base.zero() for _ in range(length)]
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            if left_index+right_index < length:
                answer[left_index+right_index] += left_value*right_value
    return tuple(answer)


relation_base = QQ if exact_mode else field
relation_ring = PolynomialRing(relation_base, names=surface_names)
monomials = tuple(
    monomial
    for degree in range(arguments.relation_degree+1)
    for monomial in relation_ring.monomials_of_degree(degree)
)
columns = []
for monomial in monomials:
    series = (relation_base.one(),)+(relation_base.zero(),)*(arguments.order-1)
    for coordinate, exponent in zip(surface_series, monomial.exponents()[0]):
        for _ in range(exponent):
            series = truncated_product(series, coordinate, arguments.order)
    columns.append(vector(field, series))
relation_matrix = Matrix(field, columns).transpose()
relation_basis = relation_matrix.right_kernel_matrix().echelon_form()
print(
    f"Q80CM24MODP|prime={'QQCM' if exact_mode else prime}|slope={arguments.slope}|"
    f"order={arguments.order}|relation_degree={arguments.relation_degree}|"
    f"monomials={len(monomials)}|rank={relation_matrix.rank()}|"
    f"kernel={relation_basis.nrows()}|status=PASS",
    flush=True,
)

common = {
    "schema": "q80-cm24-split-prime-formal-branch-v1",
    "status": "bounded_modular_formal_evidence",
    "prime": "QQCM" if exact_mode else int(prime),
    "cm_embedding": (
        {"sqrt_minus_3": str(r3), "sqrt_minus_6": str(r6)}
        if exact_mode else
        {"sqrt_minus_3": int(ZZ(r3)), "sqrt_minus_6": int(ZZ(r6))}
    ),
    "characteristic_zero_slope": arguments.slope,
    "normalization": "P=P_CM24+h",
    "order": arguments.order,
    "caveat": "finite formal reduction; not a characteristic-zero identity",
}
if arguments.jet_output and exact_mode:
    raise ValueError("exact CM jet export is intentionally disabled")
if arguments.jet_output:
    payload = dict(common)
    payload.update({
        "kind": "normalized_formal_jet",
        "active_variables": list(active_names),
        "coefficients": [[int(ZZ(value)) for value in row] for row in coefficients],
    })
    Path(arguments.jet_output).write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
if arguments.relation_basis_output:
    payload = dict(common)
    payload.update({
        "kind": "canonical_centered_relation_space",
        "centered_variables": list(surface_names),
        "relation_degree": arguments.relation_degree,
        "monomial_exponents": [list(monomial.exponents()[0]) for monomial in monomials],
        "rref_basis": [
            [str(value) if exact_mode else int(ZZ(value)) for value in row]
            for row in relation_basis.rows()
        ],
        "pivot_columns": list(relation_basis.pivots()),
    })
    Path(arguments.relation_basis_output).write_text(
        json.dumps(payload, indent=2, sort_keys=True)+"\n"
    )

Jet2.precision = 3
