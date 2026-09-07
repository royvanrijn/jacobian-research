#!/usr/bin/env sage
"""Extend and rationally reconstruct the two certified q=80 GF(7) branches.

This loads the order-seven tangent-cone certificate, continues each normalized
branch to a requested formal order, and searches for unique low-degree Pade
reconstructions of every active surface/section coordinate.  It can also
compute filtered implicit-relation spaces and export the centered two-cubic,
twelve-quartic ideal.  A reconstructed branch is accepted only when all
resolved equations vanish identically in GF(7)(h); finite-jet relations and
their exported ideal are explicitly reported only as bounded evidence. An
existing ideal can also be checked against subsequently continued, withheld
jet coefficients without rerunning the relation extraction.
"""

import argparse
import json
from pathlib import Path

from sage.all import HyperellipticCurve, Matrix, PolynomialRing, PowerSeriesRing, vector


parser = argparse.ArgumentParser()
parser.add_argument("--order", type=int, default=20)
parser.add_argument("--max-degree", type=int, default=10)
parser.add_argument(
    "--relation-degree",
    type=int,
    default=0,
    help="test implicit surface relations through this total degree",
)
parser.add_argument(
    "--pair-max-degree",
    type=int,
    default=0,
    help="scan surface-coordinate plane projections through this bidegree",
)
parser.add_argument(
    "--pair-validation",
    type=int,
    default=4,
    help="minimum withheld jet coefficients for plane relations",
)
parser.add_argument(
    "--slope",
    action="append",
    type=int,
    choices=(5, 3),
    help="restrict to one tangent slope; repeat to request both",
)
parser.add_argument(
    "--relation-summary-only",
    action="store_true",
    help="report relation dimensions without printing polynomial bases",
)
parser.add_argument(
    "--minimal-quartics",
    action="store_true",
    help="extract degree-four generators modulo affine-linear cubic multiples",
)
parser.add_argument(
    "--centered-relations",
    action="store_true",
    help="use d-3,p-4,q-3,e-2 as the four relation coordinates",
)
parser.add_argument(
    "--ideal-output",
    help="write the two cubics and twelve new quartics to a JSON artifact",
)
parser.add_argument(
    "--relation-basis-output",
    help="write the highest requested centered relation basis to a JSON artifact",
)
parser.add_argument(
    "--validate-ideal",
    help=(
        "substitute the continued surface jet into an existing ideal artifact; "
        "orders beyond its fitting order are withheld validation"
    ),
)
parser.add_argument(
    "--parameter-input",
    help="reconstruct all marked coordinates in a certified surface parameter",
)
parser.add_argument(
    "--parameter-max-degree",
    type=int,
    default=40,
)
parser.add_argument(
    "--marked-parameter-output",
    help="write a fully verified marked modular family to JSON",
)
parser.add_argument(
    "--partial-parameter-output",
    help="write individually reconstructed finite-jet coordinate candidates",
)
parser.add_argument(
    "--parameter-algebraic-coordinate",
    help="fit an implicit equation between the surface parameter and one marked coordinate",
)
parser.add_argument("--parameter-algebraic-coordinate-degree", type=int, default=2)
parser.add_argument("--parameter-algebraic-max-parameter-degree", type=int, default=40)
parser.add_argument("--parameter-algebraic-validation", type=int, default=8)
parser.add_argument("--parameter-algebraic-output")
parser.add_argument(
    "--parameter-quadratic-field-coordinate", action="append", default=[]
)
parser.add_argument("--parameter-quadratic-field-max-degree", type=int, default=30)
parser.add_argument("--parameter-quadratic-field-validation", type=int, default=12)
parser.add_argument("--parameter-quadratic-field-output")
parser.add_argument(
    "--q79-evaluation",
    type=int,
    help=(
        "evaluate Q79=-3*G1-2*G2+4*G3 at this elliptic-base value and "
        "test its direct quadratic marking cover"
    ),
)
parser.add_argument("--q79-degree-bound", type=int, default=80)
parser.add_argument("--q79-algebraic-max-parameter-degree", type=int, default=70)
parser.add_argument("--q79-validation", type=int, default=12)
parser.add_argument(
    "--q79-signs",
    default="1,1,1",
    help="independent orientations of G1,G2,G3, each +1 or -1",
)
parser.add_argument("--q79-output")
parser.add_argument(
    "--jet-input",
    help="resume a previously exported exact normalized GF(7) branch jet",
)
parser.add_argument(
    "--jet-output",
    help="export the normalized GF(7) branch jet for later continuation",
)
arguments = parser.parse_args()
if arguments.order < 8:
    parser.error("--order must be at least 8")
if arguments.max_degree < 0:
    parser.error("--max-degree must be nonnegative")
try:
    q79_signs = tuple(map(int, arguments.q79_signs.split(",")))
except ValueError:
    parser.error("--q79-signs must be a comma-separated triple of +1/-1")
if len(q79_signs) != 3 or any(sign not in (-1, 1) for sign in q79_signs):
    parser.error("--q79-signs must be a comma-separated triple of +1/-1")

load("elkies-k3/scripts/verify_q80_rank19_deformation_gf7.sage")

jet_artifact = None
if arguments.jet_input:
    jet_artifact = json.loads(Path(arguments.jet_input).read_text())
    if jet_artifact.get("schema") != "q80-cm24-normalized-formal-branch-jet-v1":
        raise ValueError("unexpected q80 branch-jet artifact schema")
    if jet_artifact.get("field") != "GF(7)":
        raise ValueError("branch-jet artifact is not over GF(7)")
    if tuple(jet_artifact.get("active_variables", ())) != active_names:
        raise ValueError("branch-jet artifact active variables do not match")


def lift_branch(slope, order_bound):
    first_order = slope*tangent[0] + tangent[1]
    assert first_order[p_column] == 1
    if jet_artifact is not None and jet_artifact["slope_mod_7"] == f"{slope}:1":
        coefficients = [
            vector(field, row) for row in jet_artifact["coefficients"]
        ]
        if len(coefficients) > order_bound:
            raise ValueError("branch-jet artifact order exceeds requested order")
        assert coefficients[0] == seed_active
        assert coefficients[1] == vector(field, first_order)
        assert coefficients[1][p_column] == 1
        assert all(not row[p_column] for row in coefficients[2:])
        assert not any(residual_coefficient(coefficients, len(coefficients)-1))
        print(
            f"Q80RANK19EXTEND|stage=jet_load|line={slope}:1|"
            f"order={len(coefficients)}|status=PASS",
            flush=True,
        )
    else:
        coefficients = [seed_active, vector(field, first_order)]
    for order in range(len(coefficients), order_bound):
        right_hand_side = -residual_coefficient(coefficients, order)
        assert not any(left_kernel*right_hand_side)
        particular = resolved_jacobian.solve_right(right_hand_side)
        # Since tangent[0] has p-coordinate zero and tangent[1] has
        # p-coordinate one, the normalization p=4+h fixes coefficient1.
        # At order >= 2 the obstruction one order later is affine in the one
        # remaining coefficient0: a nonlinear occurrence would have order at
        # least 2*order > order+1.  Two evaluations therefore replace the old
        # exhaustive seven-value scan.
        assert not tangent[0, p_column] and tangent[1, p_column] == 1
        coefficient1 = -particular[p_column]

        def normalized_trial(coefficient0):
            correction = (
                particular
                + coefficient0*tangent[0]
                + coefficient1*tangent[1]
            )
            assert not correction[p_column]
            return coefficients + [vector(field, correction)]

        trial0 = normalized_trial(field.zero())
        trial1 = normalized_trial(field.one())
        obstruction0 = left_kernel*residual_coefficient(trial0, order+1)
        obstruction1 = left_kernel*residual_coefficient(trial1, order+1)
        obstruction_step = obstruction1-obstruction0
        forced_values = {
            -constant/linear
            for constant, linear in zip(obstruction0, obstruction_step)
            if linear
        }
        assert all(
            linear or not constant
            for constant, linear in zip(obstruction0, obstruction_step)
        ), (slope, order, "inconsistent affine obstruction")
        assert len(forced_values) == 1, (slope, order, forced_values)
        coefficient0 = forced_values.pop()
        selected = normalized_trial(coefficient0)
        assert not any(
            left_kernel*residual_coefficient(selected, order+1)
        )

        # Audit the optimized selector against the previous exhaustive rule
        # at the initial orders, where the check is essentially free.
        if order <= 4:
            exhaustive = tuple(
                value for value in field
                if not any(
                    left_kernel*residual_coefficient(
                        normalized_trial(value), order+1
                    )
                )
            )
            assert exhaustive == (coefficient0,)
        coefficients = selected
        if (order + 1) % 5 == 0 or order + 1 == order_bound:
            print(
                f"Q80RANK19EXTEND|stage=lift|line={slope}:1"
                f"|order={order + 1}|unique=1",
                flush=True,
            )
    return coefficients


function_polynomials = PolynomialRing(field, "h")
h_parameter = function_polynomials.gen()
function_field = function_polynomials.fraction_field()


def pade_candidate(coefficients, numerator_degree, denominator_degree):
    """Return the unique normalized Pade candidate, if it exists."""
    count = len(coefficients)
    if denominator_degree == 0:
        if any(coefficients[index] for index in range(numerator_degree+1, count)):
            return None
        numerator = sum(
            function_polynomials(coefficients[index])*h_parameter**index
            for index in range(min(numerator_degree+1, count))
        )
        return function_field(numerator)

    row_indices = tuple(range(numerator_degree+1, count))
    if len(row_indices) < denominator_degree:
        return None
    coefficient_matrix = Matrix(
        field,
        [
            [
                coefficients[index-j] if index >= j else field.zero()
                for j in range(1, denominator_degree+1)
            ]
            for index in row_indices
        ],
    )
    target = vector(field, [-coefficients[index] for index in row_indices])
    if coefficient_matrix.rank() != denominator_degree:
        return None
    if target not in coefficient_matrix.column_space():
        return None
    denominator_tail = coefficient_matrix.solve_right(target)
    denominator_coefficients = [field.one()] + list(denominator_tail)
    numerator_coefficients = []
    for index in range(numerator_degree+1):
        numerator_coefficients.append(sum(
            denominator_coefficients[j]*coefficients[index-j]
            for j in range(min(denominator_degree, index)+1)
        ))
    for index in range(numerator_degree+1, count):
        assert not sum(
            denominator_coefficients[j]*coefficients[index-j]
            for j in range(min(denominator_degree, index)+1)
        )
    numerator = sum(
        function_polynomials(value)*h_parameter**index
        for index, value in enumerate(numerator_coefficients)
    )
    denominator = sum(
        function_polynomials(value)*h_parameter**index
        for index, value in enumerate(denominator_coefficients)
    )
    assert denominator[0] == 1
    return function_field(numerator/denominator)


def minimal_pade(coefficients, degree_bound):
    candidates = []
    for total_degree in range(2*degree_bound+1):
        for denominator_degree in range(degree_bound+1):
            numerator_degree = total_degree-denominator_degree
            if not 0 <= numerator_degree <= degree_bound:
                continue
            candidate = pade_candidate(
                coefficients, numerator_degree, denominator_degree
            )
            if candidate is not None:
                candidates.append((
                    numerator_degree,
                    denominator_degree,
                    candidate,
                ))
        if candidates:
            unique_values = {candidate for _, _, candidate in candidates}
            if len(unique_values) == 1:
                return min(candidates, key=lambda row: (max(row[:2]), row[:2]))
            return None
    return None


def bounded_pade_fast(coefficients, degree_bound):
    """Recover a rational series with numerator/denominator degree at most bound.

    Unlike ``minimal_pade``, this uses the known numerator bound to scan only
    denominator degrees.  This cuts the marked-family reconstruction from a
    quadratic grid of large column-space computations to at most ``bound``
    small consistency tests per coordinate.
    """
    count = len(coefficients)
    if count <= 2*degree_bound:
        raise ValueError("fast bounded Pade needs more than twice the degree bound")
    row_indices = tuple(range(degree_bound+1, count))
    for denominator_degree in range(degree_bound+1):
        if denominator_degree == 0:
            denominator_coefficients = [field.one()]
        else:
            coefficient_matrix = Matrix(
                field,
                [
                    [coefficients[index-j] for j in range(1, denominator_degree+1)]
                    for index in row_indices
                ],
            )
            target = vector(field, [-coefficients[index] for index in row_indices])
            if coefficient_matrix.rank() != denominator_degree:
                continue
            if coefficient_matrix.augment(target.column()).rank() != denominator_degree:
                continue
            denominator_coefficients = [field.one()] + list(
                coefficient_matrix.solve_right(target)
            )
        numerator_coefficients = [
            sum(
                denominator_coefficients[j]*coefficients[index-j]
                for j in range(min(denominator_degree, index)+1)
            )
            for index in range(count)
        ]
        if any(numerator_coefficients[degree_bound+1:]):
            continue
        numerator = sum(
            function_polynomials(value)*h_parameter**index
            for index, value in enumerate(numerator_coefficients[:degree_bound+1])
        )
        denominator = sum(
            function_polynomials(value)*h_parameter**index
            for index, value in enumerate(denominator_coefficients)
        )
        candidate = function_field(numerator/denominator)
        return (
            candidate.numerator().degree(),
            candidate.denominator().degree(),
            candidate,
        )
    return None


def tail_pade_fast(coefficients, degree_bound, verification_margin=4):
    """Recover asymmetric Pade degrees from the tail recurrence.

    A candidate of reduced degrees ``(n,d)`` is accepted from the jet only
    when ``n+d <= len(coefficients)-verification_margin``.  For each possible
    denominator degree, the final ``d+verification_margin-1`` recurrence
    equations determine and independently test its normalized denominator.
    The eventual global substitution remains the proof of exactness.
    """
    count = len(coefficients)
    candidates = []
    maximum_denominator = min(degree_bound, count-verification_margin)
    for denominator_degree in range(maximum_denominator+1):
        if denominator_degree == 0:
            denominator_coefficients = [field.one()]
        else:
            first_row = max(0, count-denominator_degree-verification_margin+1)
            row_indices = tuple(range(first_row, count))
            coefficient_matrix = Matrix(
                field,
                [
                    [
                        coefficients[index-j] if index >= j else field.zero()
                        for j in range(1, denominator_degree+1)
                    ]
                    for index in row_indices
                ],
            )
            target = vector(field, [-coefficients[index] for index in row_indices])
            if coefficient_matrix.rank() != denominator_degree:
                continue
            if coefficient_matrix.augment(target.column()).rank() != denominator_degree:
                continue
            denominator_coefficients = [field.one()] + list(
                coefficient_matrix.solve_right(target)
            )
        numerator_coefficients = [
            sum(
                denominator_coefficients[j]*coefficients[index-j]
                for j in range(min(denominator_degree, index)+1)
            )
            for index in range(count)
        ]
        nonzero_indices = [
            index for index, value in enumerate(numerator_coefficients) if value
        ]
        numerator_degree = max(nonzero_indices, default=-1)
        if numerator_degree > degree_bound:
            continue
        if numerator_degree+denominator_degree > count-verification_margin:
            continue
        numerator = sum(
            function_polynomials(value)*h_parameter**index
            for index, value in enumerate(numerator_coefficients[:numerator_degree+1])
        )
        denominator = sum(
            function_polynomials(value)*h_parameter**index
            for index, value in enumerate(denominator_coefficients)
        )
        candidate = function_field(numerator/denominator)
        candidates.append((
            candidate.numerator().degree(),
            candidate.denominator().degree(),
            candidate,
        ))
    if not candidates:
        return None
    unique_values = {candidate for _, _, candidate in candidates}
    if len(unique_values) != 1:
        return None
    return min(candidates, key=lambda row: (sum(row[:2]), max(row[:2]), row[:2]))


def reconstruct_branch(coefficients, slope):
    coordinate_series = tuple(
        tuple(row[column] for row in coefficients)
        for column in range(len(active_names))
    )
    recovered = {}
    degrees = {}
    for name, series_coefficients in zip(active_names, coordinate_series):
        result = minimal_pade(series_coefficients, arguments.max_degree)
        if result is None:
            continue
        numerator_degree, denominator_degree, candidate = result
        recovered[name] = candidate
        degrees[name] = (numerator_degree, denominator_degree)
    surface_recovered = tuple(name in recovered for name in surface_names)
    print(
        f"Q80RANK19EXTEND|stage=pade|line={slope}:1"
        f"|order={arguments.order}|max_degree={arguments.max_degree}"
        f"|recovered={len(recovered)}/{len(active_names)}"
        f"|surface={surface_recovered}",
        flush=True,
    )
    for name in surface_names:
        if name in recovered:
            print(
                f"Q80RANK19EXTEND|stage=surface_candidate|line={slope}:1"
                f"|name={name}|degrees={degrees[name]}|value={recovered[name]}",
                flush=True,
            )
    if len(recovered) != len(active_names):
        missing = tuple(name for name in active_names if name not in recovered)
        print(
            f"Q80RANK19EXTEND|stage=identity|line={slope}:1|exact=0"
            f"|reason=missing_coordinates|missing={','.join(missing)}",
            flush=True,
        )
        return False

    images = []
    for column, variable in enumerate(variables):
        if column in active_lookup:
            images.append(recovered[names[column]])
        else:
            images.append(function_field(field(seed[variable])))
    exact_map = parameters.hom(images, function_field)
    exact_values = tuple(
        exact_map(numerator)*exact_map(denominator)**(-1)
        for numerator, denominator in zip(equation_numerators, equation_denominators)
    )
    exact = not any(exact_values)
    print(
        f"Q80RANK19EXTEND|stage=identity|line={slope}:1|exact={int(exact)}",
        flush=True,
    )
    if exact:
        for name in active_names:
            print(
                f"Q80RANK19EXTEND|stage=coordinate|line={slope}:1"
                f"|name={name}|degrees={degrees[name]}|value={recovered[name]}",
                flush=True,
            )
    return exact


def reconstruct_in_surface_parameter(coefficients, slope):
    if not arguments.parameter_input:
        return None
    parameter_artifact_path = Path(arguments.parameter_input)
    parameter_artifact = json.loads(parameter_artifact_path.read_text())
    if parameter_artifact.get("schema") != "q80-cm24-formal-branch-parameter-v1":
        raise ValueError("unexpected q80 parameter artifact schema")
    if parameter_artifact.get("slope_mod_7") != f"{slope}:1":
        raise ValueError("surface parameter belongs to a different formal branch")

    parameter_polynomials = PolynomialRing(field, "t")
    t = parameter_polynomials.gen()
    parameter_function_field = parameter_polynomials.fraction_field()
    parameter_series_ring = PowerSeriesRing(
        field, "t", default_prec=arguments.order
    )
    series_t = parameter_series_ring.gen()
    centered_p_function = parameter_function_field(
        parameter_artifact["functions"]["P"]["value"]
    )
    centered_p_series = (
        parameter_series_ring(centered_p_function.numerator()(series_t))
        / parameter_series_ring(centered_p_function.denominator()(series_t))
    ).add_bigoh(arguments.order)
    if centered_p_series.valuation() != 1:
        raise AssertionError("surface parameter is not local at CM24")

    coordinate_h_series = tuple(
        tuple(row[column] for row in coefficients)
        for column in range(len(active_names))
    )
    coordinate_t_series = []
    for series_coefficients in coordinate_h_series:
        value = parameter_series_ring.zero()
        power = parameter_series_ring.one()
        for coefficient in series_coefficients:
            value += coefficient*power
            power = (power*centered_p_series).add_bigoh(arguments.order)
        coordinate_t_series.append(value.add_bigoh(arguments.order))

    # Test the geometrically marked level-79 direction directly.  In the
    # generic optimal q80 basis the exact lattice transport is
    #
    #     height4 = -G2,       Q79 = -3*G1 - 2*G2 + 4*G3.
    #
    # The separate G1,G2,G3 covers can be unnecessarily large: the Atkin--
    # Lehner marking predicts that x(Q79) descends to the rational surface
    # parameter while y(Q79) carries the genus-two quadratic extension.  A
    # single good elliptic-base evaluation is enough to test that field.
    if arguments.q79_evaluation is not None:
        base_value = field(arguments.q79_evaluation)
        if base_value in (field.zero(), field.one()):
            raise ValueError("--q79-evaluation must avoid the fixed fibers 0 and 1")
        coordinate_series_by_name = dict(zip(active_names, coordinate_t_series))

        def variable_series(name):
            if name in coordinate_series_by_name:
                return coordinate_series_by_name[name]
            variable = variables[names.index(name)]
            return parameter_series_ring(field(seed[variable])).add_bigoh(
                arguments.order
            )

        def evaluated(prefix, count):
            return sum(
                variable_series(f"{prefix}{index}") * base_value**index
                for index in range(count)
            )

        point1 = (evaluated("x1", 5), evaluated("y1", 7))
        point2 = (evaluated("x2", 5), evaluated("y2", 7))
        pole = base_value-variable_series("lam")
        if not pole[0]:
            raise ValueError("--q79-evaluation meets the G3 pole at the CM seed")
        point3 = (evaluated("n", 7)/pole**2, evaluated("m", 10)/pole**3)

        d_series = variable_series("d")
        p_series = variable_series("p")
        q_series = variable_series("q")
        r_series = -3*d_series**2+3-p_series-q_series
        a_series = base_value**2 * (
            -3 + p_series*base_value + q_series*base_value**2
            + r_series*base_value**3
        )

        def double_point(point):
            x_value, y_value = point
            if not y_value[0]:
                raise ArithmeticError("chosen q79 evaluation has singular doubling")
            slope_value = (3*x_value**2+a_series)/(2*y_value)
            x_answer = slope_value**2-2*x_value
            y_answer = slope_value*(x_value-x_answer)-y_value
            return x_answer, y_answer

        def add_distinct(left, right):
            x_left, y_left = left
            x_right, y_right = right
            if not (x_right-x_left)[0]:
                raise ArithmeticError("chosen q79 evaluation has colliding summands")
            slope_value = (y_right-y_left)/(x_right-x_left)
            x_answer = slope_value**2-x_left-x_right
            y_answer = slope_value*(x_left-x_answer)-y_left
            return x_answer, y_answer

        def negate_point(point):
            return point[0], -point[1]

        point1, point2, point3 = tuple(
            point if sign == 1 else negate_point(point)
            for point, sign in zip((point1, point2, point3), q79_signs)
        )

        # This addition chain is unit-denominator at the CM24 seed for T=5.
        minus_three_g1 = negate_point(add_distinct(double_point(point1), point1))
        minus_two_g2 = negate_point(double_point(point2))
        four_g3 = double_point(double_point(point3))
        q79_point = add_distinct(
            add_distinct(minus_three_g1, minus_two_g2), four_g3
        )
        q79_x_series, q79_y_series = q79_point

        q79_pade = tail_pade_fast(
            tuple(q79_x_series[index] for index in range(arguments.order)),
            arguments.q79_degree_bound,
            verification_margin=arguments.q79_validation,
        )
        q79_record = {
            "schema": "q80-cm24-q79-evaluation-cover-v1",
            "field": "GF(7)",
            "slope_mod_7": f"{slope}:1",
            "formal_order": int(arguments.order),
            "evaluation_T": int(base_value),
            "basis_orientation_signs": list(map(str, q79_signs)),
            "combination_coefficients_on_selected_sections": list(map(str, (
                -3*q79_signs[0], -2*q79_signs[1], 4*q79_signs[2]
            ))),
            "degree_bound": int(arguments.q79_degree_bound),
            "validation_orders": int(arguments.q79_validation),
        }
        if q79_pade is None:
            # If x(Q79) does not descend rationally in the selected marking,
            # recover its smallest quadratic field directly.  The discriminant
            # of the minimal relation determines the marked double cover and
            # can be compared with X(6,79) without reconstructing every section
            # coordinate separately.
            q79_relation_ring = PolynomialRing(field, names=("t", "x79"))
            q79_t, q79_x_variable = q79_relation_ring.gens()
            q79_x_powers = (
                parameter_series_ring.one(),
                q79_x_series,
                (q79_x_series**2).add_bigoh(arguments.order),
            )
            fitting_order = arguments.order-arguments.q79_validation
            q79_relation = None
            q79_parameter_degree = None
            for parameter_degree in range(
                arguments.q79_algebraic_max_parameter_degree+1
            ):
                monomials = tuple(
                    q79_t**parameter_exponent*q79_x_variable**x_exponent
                    for x_exponent in range(3)
                    for parameter_exponent in range(parameter_degree+1)
                )
                if len(monomials) >= fitting_order:
                    continue
                columns = []
                for monomial in monomials:
                    t_exponent, x_exponent = monomial.exponents()[0]
                    value = (
                        series_t**t_exponent*q79_x_powers[x_exponent]
                    ).add_bigoh(arguments.order)
                    columns.append(vector(
                        field, [value[index] for index in range(fitting_order)]
                    ))
                kernel = Matrix(field, columns).transpose().right_kernel_matrix()
                if kernel.nrows() != 1:
                    continue
                candidate = sum(
                    coefficient*monomial
                    for coefficient, monomial in zip(kernel.row(0), monomials)
                )
                if candidate.degree(q79_x_variable) != 2:
                    continue
                residual = parameter_series_ring.zero()
                for exponents, coefficient in candidate.dict().items():
                    residual += (
                        coefficient*series_t**exponents[0]
                        *q79_x_powers[exponents[1]]
                    )
                if residual.add_bigoh(arguments.order):
                    continue
                q79_relation = candidate
                q79_parameter_degree = parameter_degree
                break

            quadratic_branch_polynomial = None
            quadratic_genus = None
            quadratic_igusa = None
            quadratic_known_igusa = None
            quadratic_source_match = False
            if q79_relation is not None:
                coefficient_by_x_degree = tuple(
                    sum(
                        coefficient*q79_t**exponents[0]
                        for exponents, coefficient in q79_relation.dict().items()
                        if exponents[1] == x_degree
                    )
                    for x_degree in range(3)
                )
                coefficient_c, coefficient_b, coefficient_a = (
                    coefficient_by_x_degree
                )
                q79_discriminant = (
                    coefficient_b**2-4*coefficient_a*coefficient_c
                )
                quadratic_branch_polynomial = parameter_polynomials.one()
                for factor_polynomial, exponent in q79_discriminant.factor():
                    if exponent % 2:
                        quadratic_branch_polynomial *= (
                            parameter_polynomials(factor_polynomial).monic()
                        )
                quadratic_branch_degree = (
                    quadratic_branch_polynomial.degree()
                    + int(q79_discriminant.degree() % 2)
                )
                quadratic_genus = (
                    None if quadratic_branch_degree == 0
                    else (quadratic_branch_degree-2)//2
                )
                known_source = 2*parameter_polynomials.gen()**6 \
                    + 2*parameter_polynomials.gen()**4 \
                    + 4*parameter_polynomials.gen()**2 + 1
                if quadratic_genus == 2:
                    quadratic_igusa = tuple(
                        HyperellipticCurve(quadratic_branch_polynomial)
                        .absolute_igusa_invariants_kohel()
                    )
                    quadratic_known_igusa = tuple(
                        HyperellipticCurve(known_source)
                        .absolute_igusa_invariants_kohel()
                    )
                    quadratic_source_match = (
                        quadratic_igusa == quadratic_known_igusa
                    )

            q79_record.update({
                "status": (
                    "BOUNDED_QUADRATIC_X_RECONSTRUCTION"
                    if q79_relation is not None
                    else "NO_BOUNDED_DEGREE_AT_MOST_TWO_X_RECONSTRUCTION"
                ),
                "quadratic_parameter_degree": q79_parameter_degree,
                "quadratic_relation": (
                    None if q79_relation is None else str(q79_relation)
                ),
                "quadratic_squarefree_model": (
                    None if quadratic_branch_polynomial is None
                    else str(quadratic_branch_polynomial)
                ),
                "quadratic_genus": (
                    None if quadratic_genus is None else int(quadratic_genus)
                ),
                "quadratic_absolute_igusa": (
                    None if quadratic_igusa is None
                    else list(map(str, quadratic_igusa))
                ),
                "known_source_absolute_igusa": (
                    None if quadratic_known_igusa is None
                    else list(map(str, quadratic_known_igusa))
                ),
                "known_source_match": quadratic_source_match,
                "caveat": (
                    "Bounded finite-jet relation over GF(7); a global marked-"
                    "system identity and characteristic-zero lift remain pending."
                ),
            })
            print(
                f"Q80Q79COVER|line={slope}:1|T={base_value}|"
                f"signs={q79_signs}|"
                f"order={arguments.order}|degree_bound={arguments.q79_degree_bound}|"
                f"rational_x=0|quadratic_parameter_degree={q79_parameter_degree}|"
                f"quadratic_genus={quadratic_genus}|"
                f"known_source_match={int(quadratic_source_match)}|"
                "status=" + (
                    "BOUNDED_QUADRATIC_X_RECONSTRUCTION"
                    if q79_relation is not None
                    else "NO_BOUNDED_DEGREE_AT_MOST_TWO_X_RECONSTRUCTION"
                ),
                flush=True,
            )
            if q79_relation is not None:
                print(
                    f"Q80Q79COVER|quadratic_relation={q79_relation}|"
                    f"squarefree_model=w^2-({quadratic_branch_polynomial})|"
                    f"absolute_igusa={quadratic_igusa}|"
                    f"known_source_igusa={quadratic_known_igusa}",
                    flush=True,
                )
        else:
            numerator_degree, denominator_degree, q79_x_h = q79_pade
            q79_x = parameter_function_field(
                parameter_polynomials(q79_x_h.numerator().list())
                / parameter_polynomials(q79_x_h.denominator().list())
            )
            q79_x_expansion = (
                parameter_series_ring(q79_x.numerator()(series_t))
                / parameter_series_ring(q79_x.denominator()(series_t))
            ).add_bigoh(arguments.order)
            assert q79_x_expansion == q79_x_series

            centers = {"d": 3, "p": 4, "q": 3, "e": 2}
            artifact_names = {"d": "D", "p": "P", "q": "Q", "e": "E"}
            raw_surface_functions = {
                raw_name: parameter_function_field(field(center))
                + parameter_function_field(
                    parameter_artifact["functions"][artifact_names[raw_name]]["value"]
                )
                for raw_name, center in centers.items()
            }
            surface_images = [
                raw_surface_functions[name]
                if name in raw_surface_functions
                else parameter_function_field(field(seed[variable]))
                for variable, name in zip(variables, names)
            ]
            exact_surface_map = parameters.hom(
                surface_images, parameter_function_field
            )

            def map_surface_value(value):
                return (
                    exact_surface_map(value.numerator())
                    / exact_surface_map(value.denominator())
                )

            a_value = map_surface_value(A(T=base_value))
            b_value = map_surface_value(B(T=base_value))
            radicand = q79_x**3+a_value*q79_x+b_value
            assert radicand

            def odd_factors(polynomial):
                return tuple(
                    (factor_polynomial, int(exponent))
                    for factor_polynomial, exponent in polynomial.factor()
                    if exponent % 2
                )

            odd_numerator = odd_factors(radicand.numerator())
            odd_denominator = odd_factors(radicand.denominator())
            branch_polynomial = parameter_polynomials.one()
            for factor_polynomial, _ in odd_numerator+odd_denominator:
                branch_polynomial *= factor_polynomial.monic()
            finite_branch_degree = branch_polynomial.degree()
            infinity_valuation = (
                radicand.denominator().degree()-radicand.numerator().degree()
            )
            infinity_branched = bool(infinity_valuation % 2)
            branch_degree = finite_branch_degree+int(infinity_branched)
            genus = None if branch_degree == 0 else (branch_degree-2)//2
            known_source = 2*parameter_polynomials.gen()**6 \
                + 2*parameter_polynomials.gen()**4 \
                + 4*parameter_polynomials.gen()**2 + 1
            cover_igusa = None
            known_igusa = None
            known_source_match = False
            if genus == 2:
                cover_igusa = tuple(
                    HyperellipticCurve(branch_polynomial)
                    .absolute_igusa_invariants_kohel()
                )
                known_igusa = tuple(
                    HyperellipticCurve(known_source)
                    .absolute_igusa_invariants_kohel()
                )
                known_source_match = cover_igusa == known_igusa

            q79_record.update({
                "status": "BOUNDED_RATIONAL_X_RECONSTRUCTION",
                "x_degrees": [int(numerator_degree), int(denominator_degree)],
                "x": str(q79_x),
                "y_square": str(radicand),
                "odd_numerator": [
                    [str(factor), exponent] for factor, exponent in odd_numerator
                ],
                "odd_denominator": [
                    [str(factor), exponent] for factor, exponent in odd_denominator
                ],
                "squarefree_model": str(branch_polynomial),
                "infinity_valuation": int(infinity_valuation),
                "infinity_branched": infinity_branched,
                "branch_degree": int(branch_degree),
                "genus": None if genus is None else int(genus),
                "absolute_igusa": (
                    None if cover_igusa is None else list(map(str, cover_igusa))
                ),
                "known_source_absolute_igusa": (
                    None if known_igusa is None else list(map(str, known_igusa))
                ),
                "known_source_match": known_source_match,
                "caveat": (
                    "The x-coordinate is reconstructed from a finite GF(7) jet "
                    "with withheld coefficients. This is not yet a global marked-"
                    "system identity or a characteristic-zero lift."
                ),
            })
            print(
                f"Q80Q79COVER|line={slope}:1|T={base_value}|"
                f"signs={q79_signs}|"
                f"order={arguments.order}|x_degrees={numerator_degree}/{denominator_degree}|"
                f"branch_degree={branch_degree}|genus={genus}|"
                f"known_source_match={int(known_source_match)}|"
                "status=BOUNDED_RATIONAL_X_RECONSTRUCTION",
                flush=True,
            )
            print(
                f"Q80Q79COVER|squarefree_model=w^2-({branch_polynomial})|"
                f"absolute_igusa={cover_igusa}|known_source_igusa={known_igusa}",
                flush=True,
            )
        if arguments.q79_output:
            q79_path = Path(arguments.q79_output)
            q79_path.parent.mkdir(parents=True, exist_ok=True)
            q79_path.write_text(json.dumps(q79_record, indent=2, sort_keys=True)+"\n")
            print(
                f"Q80Q79COVER|output={q79_path}|status=WRITTEN",
                flush=True,
            )

    if arguments.parameter_algebraic_coordinate:
        coordinate_name = arguments.parameter_algebraic_coordinate
        if coordinate_name not in active_names:
            raise ValueError("unknown active coordinate for algebraic reconstruction")
        coordinate_series = coordinate_t_series[active_names.index(coordinate_name)]
        coordinate_degree = arguments.parameter_algebraic_coordinate_degree
        validation = arguments.parameter_algebraic_validation
        if coordinate_degree < 1 or validation < 1:
            raise ValueError("algebraic coordinate degree and validation must be positive")
        relation_ring = PolynomialRing(field, names=("t", coordinate_name))
        parameter_variable, coordinate_variable = relation_ring.gens()
        relation = None
        selected_parameter_degree = None
        for parameter_degree in range(
            arguments.parameter_algebraic_max_parameter_degree+1
        ):
            monomials = tuple(
                parameter_variable**parameter_exponent
                *coordinate_variable**coordinate_exponent
                for coordinate_exponent in range(coordinate_degree+1)
                for parameter_exponent in range(parameter_degree+1)
            )
            fitting_order = arguments.order-validation
            if len(monomials) >= fitting_order:
                continue
            columns = []
            coordinate_power = parameter_series_ring.one()
            coordinate_powers = [coordinate_power]
            for _ in range(coordinate_degree):
                coordinate_power = (
                    coordinate_power*coordinate_series
                ).add_bigoh(arguments.order)
                coordinate_powers.append(coordinate_power)
            for monomial in monomials:
                parameter_exponent, coordinate_exponent = monomial.exponents()[0]
                series = (
                    series_t**parameter_exponent
                    *coordinate_powers[coordinate_exponent]
                ).add_bigoh(arguments.order)
                columns.append(
                    vector(field, [series[index] for index in range(fitting_order)])
                )
            kernel = Matrix(field, columns).transpose().right_kernel_matrix()
            if kernel.nrows() != 1:
                continue
            candidate = sum(
                coefficient*monomial
                for coefficient, monomial in zip(kernel.row(0), monomials)
            )
            if candidate.degree(coordinate_variable) != coordinate_degree:
                continue
            residual = parameter_series_ring.zero()
            for exponents, coefficient in candidate.dict().items():
                residual += (
                    coefficient*series_t**exponents[0]
                    *coordinate_powers[exponents[1]]
                )
            residual = residual.add_bigoh(arguments.order)
            if residual:
                continue
            relation = candidate
            selected_parameter_degree = parameter_degree
            break
        print(
            f"Q80RANK19PARAM|stage=algebraic_coordinate|line={slope}:1|"
            f"coordinate={coordinate_name}|coordinate_degree={coordinate_degree}|"
            f"parameter_degree={selected_parameter_degree}|validation={validation}|"
            f"relation={relation}|status={'PASS' if relation is not None else 'NO_RELATION'}",
            flush=True,
        )
        if relation is not None and arguments.parameter_algebraic_output:
            algebraic_payload = {
                "schema": "q80-cm24-marked-coordinate-cover-v1",
                "status": "bounded_mod7_formal_relation_with_withheld_validation",
                "field": "GF(7)",
                "slope_mod_7": f"{slope}:1",
                "source_parameter": str(parameter_artifact_path),
                "formal_order": int(arguments.order),
                "validation_orders": int(validation),
                "parameter_variable": "t",
                "coordinate": coordinate_name,
                "coordinate_degree": int(coordinate_degree),
                "parameter_degree": int(selected_parameter_degree),
                "relation": str(relation),
                "caveat": (
                    "Finite formal relation over GF(7); a global marked-system "
                    "identity and characteristic-zero lift remain pending."
                ),
            }
            algebraic_path = Path(arguments.parameter_algebraic_output)
            algebraic_path.parent.mkdir(parents=True, exist_ok=True)
            algebraic_path.write_text(
                json.dumps(algebraic_payload, indent=2, sort_keys=True)+"\n"
            )
            print(
                f"Q80RANK19PARAM|stage=algebraic_artifact|output={algebraic_path}|"
                "status=WRITTEN",
                flush=True,
            )

    if arguments.parameter_quadratic_field_coordinate:
        if "lam" not in active_names:
            raise AssertionError("lambda coordinate disappeared")
        lam_series = coordinate_t_series[active_names.index("lam")]
        validation = arguments.parameter_quadratic_field_validation
        fitting_order = arguments.order-validation
        field_relations = {}
        relation_ring = PolynomialRing(field, names=("t", "lam", "z"))
        t_variable, lam_variable, z_variable = relation_ring.gens()
        for coordinate_name in arguments.parameter_quadratic_field_coordinate:
            if coordinate_name not in active_names or coordinate_name == "lam":
                raise ValueError(f"invalid quadratic-field coordinate {coordinate_name}")
            z_series = coordinate_t_series[active_names.index(coordinate_name)]
            selected = None
            selected_degree = None
            for parameter_degree in range(
                arguments.parameter_quadratic_field_max_degree+1
            ):
                monomials = tuple(
                    t_variable**parameter_exponent
                    *lam_variable**lam_exponent*z_variable**z_exponent
                    for z_exponent in range(2)
                    for lam_exponent in range(2)
                    for parameter_exponent in range(parameter_degree+1)
                )
                if len(monomials) >= fitting_order:
                    continue
                columns = []
                for monomial in monomials:
                    t_exponent, lam_exponent, z_exponent = monomial.exponents()[0]
                    value = series_t**t_exponent
                    if lam_exponent:
                        value *= lam_series
                    if z_exponent:
                        value *= z_series
                    value = value.add_bigoh(arguments.order)
                    columns.append(
                        vector(field, [value[index] for index in range(fitting_order)])
                    )
                kernel = Matrix(field, columns).transpose().right_kernel_matrix()
                if kernel.nrows() != 1:
                    continue
                candidate = sum(
                    coefficient*monomial
                    for coefficient, monomial in zip(kernel.row(0), monomials)
                )
                if candidate.degree(z_variable) != 1:
                    continue
                residual = parameter_series_ring.zero()
                for exponents, coefficient in candidate.dict().items():
                    value = coefficient*series_t**exponents[0]
                    if exponents[1]:
                        value *= lam_series
                    if exponents[2]:
                        value *= z_series
                    residual += value
                if residual.add_bigoh(arguments.order):
                    continue
                selected = candidate
                selected_degree = parameter_degree
                break
            field_relations[coordinate_name] = {
                "parameter_degree": selected_degree,
                "relation": None if selected is None else str(selected),
            }
            print(
                f"Q80RANK19PARAM|stage=quadratic_field_coordinate|line={slope}:1|"
                f"coordinate={coordinate_name}|parameter_degree={selected_degree}|"
                f"validation={validation}|status={'PASS' if selected is not None else 'NO_RELATION'}",
                flush=True,
            )
        if arguments.parameter_quadratic_field_output:
            output_path = Path(arguments.parameter_quadratic_field_output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(json.dumps({
                "schema": "q80-cm24-p3-quadratic-field-coordinates-v1",
                "status": "bounded_mod7_formal_relations_with_withheld_validation",
                "field": "GF(7)",
                "slope_mod_7": f"{slope}:1",
                "source_parameter": str(parameter_artifact_path),
                "formal_order": int(arguments.order),
                "validation_orders": int(validation),
                "quadratic_generator": "lam",
                "coordinates": field_relations,
                "caveat": "Finite formal relations; global marked identities remain pending.",
            }, indent=2, sort_keys=True)+"\n")
            print(
                f"Q80RANK19PARAM|stage=quadratic_field_artifact|output={output_path}|"
                "status=WRITTEN",
                flush=True,
            )

    recovered = {}
    degrees = {}
    for name, series in zip(active_names, coordinate_t_series):
        result = tail_pade_fast(
            tuple(series[index] for index in range(arguments.order)),
            arguments.parameter_max_degree,
        )
        if result is None:
            continue
        numerator_degree, denominator_degree, h_candidate = result
        numerator = parameter_polynomials(h_candidate.numerator().list())
        denominator = parameter_polynomials(h_candidate.denominator().list())
        recovered[name] = parameter_function_field(numerator/denominator)
        degrees[name] = (numerator_degree, denominator_degree)

    print(
        f"Q80RANK19PARAM|stage=pade|line={slope}:1|order={arguments.order}|"
        f"max_degree={arguments.parameter_max_degree}|"
        f"recovered={len(recovered)}/{len(active_names)}",
        flush=True,
    )
    print(
        f"Q80RANK19PARAM|stage=recovered_degrees|line={slope}:1|values="
        + ",".join(
            f"{name}:{degrees[name][0]}/{degrees[name][1]}"
            for name in active_names if name in degrees
        ),
        flush=True,
    )
    for name in active_names:
        if name in recovered:
            print(
                f"Q80RANK19PARAM|stage=recovered_coordinate|line={slope}:1|"
                f"name={name}|degrees={degrees[name]}|value={recovered[name]}",
                flush=True,
            )
    if arguments.partial_parameter_output:
        partial_artifact = {
            "schema": "q80-cm24-marked-branch-partial-parameter-v1",
            "status": "bounded_mod7_formal_coordinate_candidates",
            "field": "GF(7)",
            "slope_mod_7": f"{slope}:1",
            "source_parameter": str(parameter_artifact_path),
            "formal_order": int(arguments.order),
            "coordinates": {
                name: {
                    "degrees": list(map(int, degrees[name])),
                    "value": str(recovered[name]),
                }
                for name in active_names if name in recovered
            },
            "caveat": (
                "Each entry is a unique bounded finite-jet candidate. Only a "
                "subsequent global identity check may promote it to exact."
            ),
        }
        partial_path = Path(arguments.partial_parameter_output)
        partial_path.parent.mkdir(parents=True, exist_ok=True)
        partial_path.write_text(
            json.dumps(partial_artifact, indent=2, sort_keys=True)+"\n"
        )
        print(
            f"Q80RANK19PARAM|stage=partial_artifact|output={partial_path}|"
            f"coordinates={len(recovered)}|status=WRITTEN",
            flush=True,
        )
    if len(recovered) != len(active_names):
        missing = tuple(name for name in active_names if name not in recovered)
        print(
            f"Q80RANK19PARAM|stage=identity|line={slope}:1|exact=0|"
            f"reason=missing_coordinates|missing={','.join(missing)}",
            flush=True,
        )
        return False

    images = []
    for column, variable in enumerate(variables):
        if column in active_lookup:
            images.append(recovered[names[column]])
        else:
            images.append(parameter_function_field(field(seed[variable])))
    exact_map = parameters.hom(images, parameter_function_field)
    exact_values = tuple(
        exact_map(numerator)*exact_map(denominator)**(-1)
        for numerator, denominator in zip(equation_numerators, equation_denominators)
    )
    exact = not any(exact_values)

    center_by_name = {"d": field(3), "p": field(4), "q": field(3), "e": field(2)}
    artifact_by_name = {"d": "D", "p": "P", "q": "Q", "e": "E"}
    surface_match = True
    for raw_name, artifact_name in artifact_by_name.items():
        expected = parameter_function_field(
            parameter_artifact["functions"][artifact_name]["value"]
        )
        surface_match &= recovered[raw_name]-center_by_name[raw_name] == expected
    exact &= surface_match
    print(
        f"Q80RANK19PARAM|stage=identity|line={slope}:1|exact={int(exact)}|"
        f"surface_artifact_match={int(surface_match)}",
        flush=True,
    )
    if exact:
        for name in active_names:
            print(
                f"Q80RANK19PARAM|stage=coordinate|line={slope}:1|name={name}|"
                f"degrees={degrees[name]}|value={recovered[name]}",
                flush=True,
            )
        if arguments.marked_parameter_output:
            output_artifact = {
                "schema": "q80-cm24-marked-branch-parameter-v1",
                "status": "exact_mod7_identity_on_resolved_marked_system",
                "field": "GF(7)",
                "slope_mod_7": f"{slope}:1",
                "characteristic_zero_slope": parameter_artifact[
                    "characteristic_zero_slope"
                ],
                "source_parameter": str(parameter_artifact_path),
                "formal_order": int(arguments.order),
                "coordinates": {
                    name: {
                        "degrees": list(map(int, degrees[name])),
                        "value": str(recovered[name]),
                    }
                    for name in active_names
                },
                "verification": {
                    "all_resolved_equations": True,
                    "surface_parameter_match": True,
                },
                "caveat": (
                    "This is an exact identity over GF(7)(t) for the resolved marked "
                    "system; a characteristic-zero lift remains open."
                ),
            }
            output_path = Path(arguments.marked_parameter_output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(
                json.dumps(output_artifact, indent=2, sort_keys=True)+"\n"
            )
            print(
                f"Q80RANK19PARAM|stage=artifact|output={output_path}|status=WRITTEN",
                flush=True,
            )
    return exact


surface_relation_ring = PolynomialRing(field, names=("d", "p", "q", "e"))


def truncated_product(left, right):
    result = [field.zero() for _ in range(arguments.order)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            if i+j < arguments.order:
                result[i+j] += left_value*right_value
    return tuple(result)


def report_surface_relations(coefficients, slope):
    if not arguments.relation_degree:
        return
    surface_series = tuple(
        tuple(row[column] for row in coefficients) for column in surface_columns
    )
    if arguments.centered_relations:
        surface_series = tuple(
            (field.zero(),)+series[1:] for series in surface_series
        )
    relation_coordinates = "centered_at_CM24" if arguments.centered_relations else "raw"
    relation_data = {}
    for degree in range(1, arguments.relation_degree+1):
        monomials = tuple(
            monomial
            for current_degree in range(degree+1)
            for monomial in surface_relation_ring.monomials_of_degree(current_degree)
        )
        columns = []
        for monomial in monomials:
            series = (field.one(),)+(field.zero(),)*(arguments.order-1)
            for coordinate, exponent in zip(surface_series, monomial.exponents()[0]):
                for _ in range(exponent):
                    series = truncated_product(series, coordinate)
            columns.append(vector(field, series))
        kernel = Matrix(field, columns).transpose().right_kernel_matrix()
        relation_data[degree] = (monomials, kernel)
        determined = arguments.order > len(monomials)
        relations = tuple(
            sum(coefficient*monomial for coefficient, monomial in zip(row, monomials))
            for row in kernel.rows()
        ) if determined and not arguments.relation_summary_only else ()
        print(
            f"Q80RANK19RELATION|line={slope}:1|degree={degree}|"
            f"order={arguments.order}|monomials={len(monomials)}|"
            f"determined={int(determined)}|kernel={kernel.nrows()}|"
            f"quotient_dimension={len(monomials)-kernel.nrows()}|"
            f"coordinates={relation_coordinates}|"
            f"relations={relations}",
            flush=True,
        )
    if arguments.relation_basis_output:
        if not arguments.centered_relations:
            raise ValueError("--relation-basis-output requires --centered-relations")
        degree = arguments.relation_degree
        monomials, kernel = relation_data[degree]
        if arguments.order <= len(monomials) or not kernel.nrows():
            raise ValueError("highest requested relation space is not determined and nonzero")
        relations = tuple(
            sum(
                coefficient*monomial
                for coefficient, monomial in zip(row, monomials)
            )
            for row in kernel.rows()
        )
        affine_export_ring = PolynomialRing(field, names=("D", "P", "Q", "E"))
        affine_generators = tuple(affine_export_ring({
            exponents: coefficient
            for exponents, coefficient in relation.dict().items()
        }) for relation in relations)
        artifact = {
            "schema": "q80-cm24-formal-branch-ideal-v1",
            "status": "bounded_mod7_formal_evidence",
            "field": "GF(7)",
            "slope_mod_7": f"{slope}:1",
            "characteristic_zero_slope": "8/87" if slope == 5 else "1/12",
            "formal_order": arguments.order,
            "raw_center": {"d": 3, "p": 4, "q": 3, "e": 2},
            "affine_variables": ["D", "P", "Q", "E"],
            "affine_coordinate_definition": {
                "D": "d-3", "P": "p-4", "Q": "q-3", "E": "e-2"
            },
            "affine_generators": list(map(str, affine_generators)),
            "generator_degrees": [degree]*len(affine_generators),
            "relation_space": {
                "degree": degree,
                "monomials": len(monomials),
                "kernel_dimension": kernel.nrows(),
                "surplus_jet_equations": arguments.order-len(monomials),
            },
            "caveat": (
                "The generators annihilate the selected finite formal jet. "
                "Continued-jet validation and global algebraization remain separate."
            ),
        }
        output_path = Path(arguments.relation_basis_output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(artifact, indent=2, sort_keys=True, default=int)+"\n"
        )
        print(
            f"Q80RANK19IDEAL|line={slope}:1|degree={degree}|"
            f"generators={len(affine_generators)}|output={output_path}|status=WRITTEN",
            flush=True,
        )
    if arguments.minimal_quartics:
        if arguments.relation_degree < 4:
            raise ValueError("--minimal-quartics requires --relation-degree at least 4")
        cubic_monomials, cubic_kernel = relation_data[3]
        quartic_monomials, quartic_kernel = relation_data[4]
        cubic_relations = tuple(
            sum(coefficient*monomial for coefficient, monomial in zip(row, cubic_monomials))
            for row in cubic_kernel.rows()
        )
        ring_variables = surface_relation_ring.gens()

        def coefficient_row(polynomial):
            return vector(
                field,
                [polynomial.monomial_coefficient(monomial) for monomial in quartic_monomials],
            )

        multiple_rows = Matrix(
            field,
            [coefficient_row(multiplier*relation)
             for relation in cubic_relations
             for multiplier in (surface_relation_ring.one(),)+ring_variables],
        )
        assert multiple_rows.rank() == 10
        multiple_echelon = multiple_rows.echelon_form()
        multiple_echelon = Matrix(
            field, [row for row in multiple_echelon.rows() if row]
        )
        multiple_pivots = multiple_echelon.pivots()
        remainders = []
        for row in quartic_kernel.rows():
            remainder = vector(field, row)
            for index, pivot in enumerate(multiple_pivots):
                remainder -= remainder[pivot]*multiple_echelon[index]
            remainders.append(remainder)
        new_echelon = Matrix(field, remainders).echelon_form()
        new_rows = tuple(row for row in new_echelon.rows() if row)
        assert len(new_rows) == 12
        new_relations = tuple(
            sum(coefficient*monomial for coefficient, monomial in zip(row, quartic_monomials))
            for row in new_rows
        )
        support_sizes = tuple(len(relation.monomials()) for relation in new_relations)
        reported_new_relations = () if arguments.relation_summary_only else new_relations
        print(
            f"Q80RANK19MINIMAL|line={slope}:1|cubic_generators=2|"
            f"cubic_linear_multiples=10|new_quartic_generators={len(new_relations)}|"
            f"coordinates={relation_coordinates}|support_sizes={support_sizes}|"
            f"relations={reported_new_relations}",
            flush=True,
        )
        if arguments.ideal_output:
            if not arguments.centered_relations:
                raise ValueError("--ideal-output requires --centered-relations")
            affine_export_ring = PolynomialRing(
                field, names=("D", "P", "Q", "E")
            )
            homogeneous_export_ring = PolynomialRing(
                field, names=("z", "D", "P", "Q", "E")
            )
            z, *homogeneous_variables = homogeneous_export_ring.gens()

            def export_affine(polynomial):
                return affine_export_ring({
                    exponents: coefficient
                    for exponents, coefficient in polynomial.dict().items()
                })

            def export_homogeneous(polynomial, total_degree):
                result = homogeneous_export_ring.zero()
                for exponents, coefficient in polynomial.dict().items():
                    monomial = z**(total_degree-sum(exponents))
                    for variable, exponent in zip(homogeneous_variables, exponents):
                        monomial *= variable**exponent
                    result += coefficient*monomial
                return result

            affine_generators = tuple(map(export_affine, cubic_relations+new_relations))
            homogeneous_generators = tuple(
                export_homogeneous(relation, total_degree)
                for relation, total_degree in zip(
                    cubic_relations+new_relations, (3, 3)+(4,)*12
                )
            )
            artifact = {
                "schema": "q80-cm24-formal-branch-ideal-v1",
                "status": "bounded_mod7_formal_evidence",
                "field": "GF(7)",
                "slope_mod_7": f"{slope}:1",
                "characteristic_zero_slope": "8/87" if slope == 5 else "1/12",
                "formal_order": arguments.order,
                "raw_center": {"d": 3, "p": 4, "q": 3, "e": 2},
                "affine_variables": ["D", "P", "Q", "E"],
                "affine_coordinate_definition": {
                    "D": "d-3", "P": "p-4", "Q": "q-3", "E": "e-2"
                },
                "affine_generators": [str(relation) for relation in affine_generators],
                "homogeneous_variables": ["z", "D", "P", "Q", "E"],
                "homogeneous_generators": [
                    str(relation) for relation in homogeneous_generators
                ],
                "generator_degrees": [3, 3]+[4]*12,
                "filtered_quotient_dimensions_through_degree_4": [1, 5, 15, 33, 48],
                "caveat": (
                    "Generators annihilate the selected formal jet through the stated "
                    "order; this artifact alone does not prove the global characteristic-"
                    "zero branch or its normalization."
                ),
            }
            output_path = Path(arguments.ideal_output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(
                json.dumps(artifact, indent=2, sort_keys=True, default=int)+"\n"
            )
            print(
                f"Q80RANK19IDEAL|line={slope}:1|output={output_path}|"
                "affine_generators=14|homogeneous_generators=14|status=WRITTEN",
                flush=True,
            )


def report_surface_pair_relations(coefficients, slope):
    if not arguments.pair_max_degree:
        return
    surface_series = tuple(
        tuple(row[column] for row in coefficients) for column in surface_columns
    )
    if arguments.centered_relations:
        surface_series = tuple(
            (field.zero(),)+series[1:] for series in surface_series
        )
    powers = []
    for coordinate in surface_series:
        coordinate_powers = [
            (field.one(),)+(field.zero(),)*(arguments.order-1)
        ]
        for _ in range(arguments.pair_max_degree):
            coordinate_powers.append(
                truncated_product(coordinate_powers[-1], coordinate)
            )
        powers.append(tuple(coordinate_powers))
    candidates = []
    for left in range(len(surface_names)):
        for right in range(left+1, len(surface_names)):
            for left_degree in range(1, arguments.pair_max_degree+1):
                for right_degree in range(1, arguments.pair_max_degree+1):
                    monomial_count = (left_degree+1)*(right_degree+1)
                    validation = arguments.order-monomial_count
                    if validation < arguments.pair_validation:
                        continue
                    columns = []
                    for left_exponent in range(left_degree+1):
                        for right_exponent in range(right_degree+1):
                            columns.append(vector(field, truncated_product(
                                powers[left][left_exponent],
                                powers[right][right_exponent],
                            )))
                    kernel_dimension = Matrix(field, columns).transpose().right_kernel().dimension()
                    if kernel_dimension:
                        candidates.append((
                            surface_names[left],
                            surface_names[right],
                            left_degree,
                            right_degree,
                            validation,
                            kernel_dimension,
                        ))
    print(
        f"Q80RANK19PAIR|line={slope}:1|order={arguments.order}|"
        f"max_degree={arguments.pair_max_degree}|"
        f"validation_min={arguments.pair_validation}|"
        f"candidates={candidates}|status=PASS",
        flush=True,
    )


def validate_surface_ideal(coefficients, slope):
    if not arguments.validate_ideal:
        return True
    artifact_path = Path(arguments.validate_ideal)
    artifact = json.loads(artifact_path.read_text())
    if artifact.get("schema") != "q80-cm24-formal-branch-ideal-v1":
        raise ValueError("unexpected q80 ideal artifact schema")
    if artifact.get("field") != "GF(7)":
        raise ValueError("ideal artifact is not over GF(7)")
    if artifact.get("slope_mod_7") != f"{slope}:1":
        raise ValueError("ideal artifact belongs to a different formal branch")

    ideal_ring = PolynomialRing(field, names=tuple(artifact["affine_variables"]))
    generators = tuple(
        ideal_ring(polynomial) for polynomial in artifact["affine_generators"]
    )
    centered_surface_series = tuple(
        (field.zero(),) + tuple(row[column] for row in coefficients[1:])
        for column in surface_columns
    )

    def evaluate_generator(polynomial):
        result = (field.zero(),) * arguments.order
        for exponents, coefficient in polynomial.dict().items():
            term = (field(coefficient),) + (field.zero(),) * (arguments.order - 1)
            for coordinate, exponent in zip(centered_surface_series, exponents):
                for _ in range(exponent):
                    term = truncated_product(term, coordinate)
            result = tuple(left + right for left, right in zip(result, term))
        return result

    residuals = tuple(evaluate_generator(generator) for generator in generators)
    failures = tuple(
        (generator_index, coefficient_index, coefficient)
        for generator_index, residual in enumerate(residuals)
        for coefficient_index, coefficient in enumerate(residual)
        if coefficient
    )
    fitted_order = int(artifact["formal_order"])
    withheld = max(arguments.order - fitted_order, 0)
    print(
        f"Q80RANK19IDEALVALIDATE|line={slope}:1|generators={len(generators)}|"
        f"artifact_order={fitted_order}|continued_order={arguments.order}|"
        f"withheld_coefficients={withheld}|failures={failures[:8]}|"
        f"status={'PASS' if not failures else 'FAIL'}",
        flush=True,
    )
    return not failures


exact_branches = 0
parameter_exact_branches = 0
selected_slopes = tuple(arguments.slope) if arguments.slope else (5, 3)
for branch_slope in selected_slopes:
    branch_coefficients = lift_branch(branch_slope, arguments.order)
    if arguments.jet_output:
        if len(selected_slopes) != 1:
            raise ValueError("--jet-output requires exactly one selected slope")
        jet_path = Path(arguments.jet_output)
        jet_path.parent.mkdir(parents=True, exist_ok=True)
        jet_path.write_text(json.dumps({
            "schema": "q80-cm24-normalized-formal-branch-jet-v1",
            "status": "exact_finite_formal_jet",
            "field": "GF(7)",
            "slope_mod_7": f"{branch_slope}:1",
            "characteristic_zero_slope": (
                "8/87" if branch_slope == 5 else "1/12"
            ),
            "normalization": "p=4+h",
            "order": len(branch_coefficients),
            "active_variables": list(active_names),
            "coefficients": [
                list(map(int, coefficient))
                for coefficient in branch_coefficients
            ],
            "caveat": (
                "This is an exact finite formal jet over GF(7), not a global "
                "algebraization or characteristic-zero identity."
            ),
        }, indent=2, sort_keys=True)+"\n")
        print(
            f"Q80RANK19EXTEND|stage=jet_output|line={branch_slope}:1|"
            f"order={len(branch_coefficients)}|output={jet_path}|status=WRITTEN",
            flush=True,
        )
    report_surface_relations(branch_coefficients, branch_slope)
    report_surface_pair_relations(branch_coefficients, branch_slope)
    if not validate_surface_ideal(branch_coefficients, branch_slope):
        raise AssertionError("candidate ideal fails continued-jet validation")
    exact_branches += int(reconstruct_branch(branch_coefficients, branch_slope))
    parameter_result = reconstruct_in_surface_parameter(
        branch_coefficients, branch_slope
    )
    parameter_exact_branches += int(bool(parameter_result))

print(
    f"Q80RANK19EXTEND|order={arguments.order}|max_degree={arguments.max_degree}"
    f"|branches={len(selected_slopes)}|exact_branches={exact_branches}|"
    f"parameter_exact_branches={parameter_exact_branches}|status=PASS",
    flush=True,
)
