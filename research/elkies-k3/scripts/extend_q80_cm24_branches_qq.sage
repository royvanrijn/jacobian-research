#!/usr/bin/env sage
"""Extend the two exact q=80 CM24 rank-19 branches over characteristic zero.

This loads the exact CM24 model and its resolved two-dimensional tangent cone,
normalizes the local parameter by ``P=9/4+h``, and recursively solves the
marked section equations over ``Q(sqrt(-3),sqrt(-6))[[h]]``.  At each order,
the next obstruction fixes the one remaining tangent-kernel ambiguity.  Pade
candidates and implicit relations are finite-jet diagnostics; only the
displayed formal coefficients and exact residual checks through the requested
order are certified here.  For the `8/87` cubic space, the script also checks
reduction against the independently continued order-145 GF(7) reference span.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, GCD, LCM, Matrix, PolynomialRing, QQ, ZZ, vector


parser = argparse.ArgumentParser()
parser.add_argument("--order", type=int, default=14)
parser.add_argument("--max-degree", type=int, default=7)
parser.add_argument("--relation-degree", type=int, default=2)
parser.add_argument("--pair-max-degree", type=int, default=5)
parser.add_argument("--pair-validation", type=int, default=3)
parser.add_argument(
    "--pair-only",
    nargs=2,
    choices=("D", "P", "Q", "E"),
    metavar=("X", "Y"),
    help="restrict implicit pair diagnostics to one ordered coordinate pair",
)
parser.add_argument(
    "--pair-total-degree",
    type=int,
    default=0,
    help="test one total-degree plane relation instead of bidegree rectangles",
)
parser.add_argument("--quiet-series", action="store_true")
parser.add_argument(
    "--surface-series-output",
    help="write the exact QQ surface series to a generated JSON artifact",
)
parser.add_argument(
    "--surface-series-only",
    action="store_true",
    help="skip relation, pair, and Pade diagnostics after lifting the branch",
)
parser.add_argument(
    "--supported-plane-ideal",
    help="use a GF(7) ideal artifact to recover a supported exact plane relation",
)
parser.add_argument(
    "--supported-plane-pair",
    nargs=2,
    default=("D", "Q"),
    metavar=("X", "Y"),
)
parser.add_argument(
    "--slope",
    action="append",
    choices=("8/87", "1/12"),
    help="restrict to one tangent slope; repeat to request both",
)
arguments = parser.parse_args()
if arguments.order < 4:
    parser.error("--order must be at least 4")
if arguments.pair_validation < 1:
    parser.error("--pair-validation must be positive")

load("elkies-k3/scripts/verify_q80_cm24_rational_model.sage")

seed_active = vector(
    composite, [cone_seed_values[column] for column in cone_active_columns]
)
surface_names = ("D", "P", "Q", "E")
surface_columns = tuple(cone_active_names.index(name) for name in surface_names)
p_column = cone_active_names.index("P")

# The same rank-37 linear system is solved at every formal order.  Pin two
# free coordinates to zero and invert one nonsingular 37-by-37 minor once,
# rather than asking Sage to repeat a number-field elimination each time.
solver_columns = tuple(cone_resolved_jacobian.pivots())
solver_column_basis = cone_resolved_jacobian.matrix_from_columns(solver_columns)
solver_rows = tuple(solver_column_basis.transpose().pivots())
solver_minor_inverse = solver_column_basis.matrix_from_rows(solver_rows).inverse()
assert len(solver_columns) == len(solver_rows) == 37


def solve_particular(right_hand_side):
    solution = vector(composite, [0]*cone_resolved_jacobian.ncols())
    pivot_values = solver_minor_inverse*vector(
        composite, [right_hand_side[index] for index in solver_rows]
    )
    for column, value in zip(solver_columns, pivot_values):
        solution[column] = value
    assert cone_resolved_jacobian*solution == right_hand_side
    return solution


class FastCM:
    """Four rational coordinates in QQ(sqrt(-3),sqrt(-6))."""

    __slots__ = ("coordinates",)

    def __init__(self, a=0, b=0, c=0, d=0):
        self.coordinates = (QQ(a), QQ(b), QQ(c), QQ(d))

    def __bool__(self):
        return any(self.coordinates)

    def __hash__(self):
        return hash(self.coordinates)

    def __eq__(self, other):
        other = fast_cm(other)
        return self.coordinates == other.coordinates

    def __repr__(self):
        return f"FastCM{self.coordinates}"

    def __add__(self, other):
        other = fast_cm(other)
        return FastCM(*(a+b for a, b in zip(self.coordinates, other.coordinates)))

    __radd__ = __add__

    def __neg__(self):
        return FastCM(*(-value for value in self.coordinates))

    def __sub__(self, other):
        return self+(-fast_cm(other))

    def __rsub__(self, other):
        return fast_cm(other)-self

    def __mul__(self, other):
        other = fast_cm(other)
        a, b, c, d = self.coordinates
        A, B, C, D = other.coordinates
        return FastCM(
            a*A-3*b*B-6*c*C+18*d*D,
            a*B+b*A-6*c*D-6*d*C,
            a*C+c*A-3*b*D-3*d*B,
            a*D+d*A+b*C+c*B,
        )

    __rmul__ = __mul__

    def conjugate_three(self):
        a, b, c, d = self.coordinates
        return FastCM(a, -b, c, -d)

    def conjugate_six(self):
        a, b, c, d = self.coordinates
        return FastCM(a, b, -c, -d)

    def conjugate_both(self):
        a, b, c, d = self.coordinates
        return FastCM(a, -b, -c, d)

    def inverse(self):
        numerator = (
            self.conjugate_three()
            *self.conjugate_six()
            *self.conjugate_both()
        )
        norm = self*numerator
        assert not any(norm.coordinates[1:]) and norm.coordinates[0]
        return FastCM(*(value/norm.coordinates[0] for value in numerator.coordinates))

    def __truediv__(self, other):
        return self*fast_cm(other).inverse()

    def __rtruediv__(self, other):
        return fast_cm(other)*self.inverse()

    def rational(self):
        assert not any(self.coordinates[1:])
        return self.coordinates[0]


basis_matrix_inverse = basis_matrix.inverse()


def fast_cm(value):
    if isinstance(value, FastCM):
        return value
    try:
        if value.parent() is composite:
            coordinates = basis_matrix_inverse*vector(QQ, vector(composite(value)))
            return FastCM(*coordinates)
    except AttributeError:
        pass
    return FastCM(QQ(value))


class FastCMRing:
    def __call__(self, value=0):
        return fast_cm(value)

    @staticmethod
    def zero():
        return FastCM()

    @staticmethod
    def one():
        return FastCM(1)


fast_cm_ring = FastCMRing()


def fast_vector(values):
    return [fast_cm(value) for value in values]


def fast_matrix(rows):
    return [[fast_cm(value) for value in row] for row in rows]


def fast_add(left, right):
    return [a+b for a, b in zip(left, right)]


def fast_subtract(left, right):
    return [a-b for a, b in zip(left, right)]


def fast_scale(scalar, values):
    scalar = fast_cm(scalar)
    return [scalar*value for value in values]


def fast_matvec(matrix_rows, values):
    return [
        sum((coefficient*value for coefficient, value in zip(row, values)), FastCM())
        for row in matrix_rows
    ]


fast_seed_active = fast_vector(seed_active)
fast_tangent = [fast_vector(row) for row in cone_tangent.rows()]
fast_solver_inverse = fast_matrix(solver_minor_inverse.rows())
fast_solver_column_basis = fast_matrix(solver_column_basis.rows())
fast_left_kernel = fast_matrix(cone_left_kernel.rows())


def fast_solve_particular(right_hand_side):
    selected = [right_hand_side[index] for index in solver_rows]
    pivot_values = fast_matvec(fast_solver_inverse, selected)
    solution = [FastCM() for _ in range(cone_resolved_jacobian.ncols())]
    for column, value in zip(solver_columns, pivot_values):
        solution[column] = value
    assert fast_matvec(fast_solver_column_basis, pivot_values) == right_hand_side
    return solution


def fast_residual_values(coefficients, precision):
    Jet2.precision = precision
    Jet2.base_ring = fast_cm_ring
    active_lookup = {
        column: index for index, column in enumerate(cone_active_columns)
    }
    inputs = []
    for column, seed_value in enumerate(cone_seed_values):
        if column in active_lookup:
            active_index = active_lookup[column]
            series = [coefficient[active_index] for coefficient in coefficients]
        else:
            series = [fast_cm(seed_value)]
        inputs.append(Jet2.from_coefficients(series))
    return evaluate_cone_system(inputs)


def fast_residual_coefficient(coefficients, order):
    return [
        value.coefficients[order]
        for value in fast_residual_values(coefficients, order+1)
    ]


def lift_branch_fast(slope):
    first_order = fast_add(
        fast_scale(slope, fast_tangent[0]), fast_tangent[1]
    )
    assert first_order[p_column] == 1
    coefficients = [fast_seed_active, first_order]
    for order in range(2, arguments.order):
        residual = fast_residual_coefficient(coefficients, order)
        assert not any(fast_matvec(fast_left_kernel, residual))
        particular = fast_solve_particular(fast_scale(-1, residual))
        correction0 = fast_subtract(
            particular,
            fast_scale(particular[p_column], fast_tangent[1]),
        )
        assert not correction0[p_column]
        trial0 = coefficients+[correction0]
        trial1 = coefficients+[fast_add(correction0, fast_tangent[0])]
        obstruction0 = fast_matvec(
            fast_left_kernel, fast_residual_coefficient(trial0, order+1)
        )
        obstruction1 = fast_matvec(
            fast_left_kernel, fast_residual_coefficient(trial1, order+1)
        )
        obstruction_delta = fast_subtract(obstruction1, obstruction0)
        alpha_values = {
            -left/right
            for left, right in zip(obstruction0, obstruction_delta)
            if right
        }
        assert len(alpha_values) == 1
        alpha = alpha_values.pop()
        assert not any(
            left+alpha*right
            for left, right in zip(obstruction0, obstruction_delta)
        )
        coefficients.append(
            fast_add(correction0, fast_scale(alpha, fast_tangent[0]))
        )
    assert not any(
        coefficient
        for residual in fast_residual_values(coefficients, arguments.order)
        for coefficient in residual.coefficients
    )
    Jet2.base_ring = composite
    return coefficients


def residual_values(coefficients, precision):
    Jet2.precision = precision
    active_lookup = {
        column: index for index, column in enumerate(cone_active_columns)
    }
    inputs = []
    for column, seed_value in enumerate(cone_seed_values):
        if column in active_lookup:
            active_index = active_lookup[column]
            series = [coefficient[active_index] for coefficient in coefficients]
        else:
            series = [seed_value]
        inputs.append(Jet2.from_coefficients(series))
    return evaluate_cone_system(inputs)


def residual_coefficient(coefficients, order):
    return vector(
        composite,
        [value.coefficients[order] for value in residual_values(coefficients, order+1)],
    )


def lift_branch(slope):
    first_order = slope*cone_tangent[0]+cone_tangent[1]
    assert first_order[p_column] == 1
    coefficients = [seed_active, vector(composite, first_order)]
    for order in range(2, arguments.order):
        residual = residual_coefficient(coefficients, order)
        assert not any(cone_left_kernel*residual)
        particular = solve_particular(-residual)
        # P=9/4+h removes reparametrization.  Since v0.P=0 and v1.P=1,
        # the v1 coefficient is fixed immediately and only v0 remains.
        correction0 = particular-particular[p_column]*cone_tangent[1]
        assert correction0[p_column] == 0
        trial0 = coefficients+[vector(composite, correction0)]
        trial1 = coefficients+[vector(composite, correction0+cone_tangent[0])]
        obstruction0 = cone_left_kernel*residual_coefficient(trial0, order+1)
        obstruction_delta = (
            cone_left_kernel*residual_coefficient(trial1, order+1)-obstruction0
        )
        alpha_values = {
            -left/right
            for left, right in zip(obstruction0, obstruction_delta)
            if right
        }
        assert len(alpha_values) == 1
        alpha = alpha_values.pop()
        assert not any(
            left+alpha*right
            for left, right in zip(obstruction0, obstruction_delta)
        )
        coefficients.append(
            vector(composite, correction0+alpha*cone_tangent[0])
        )
    full_residuals = residual_values(coefficients, arguments.order)
    assert not any(
        coefficient
        for residual in full_residuals
        for coefficient in residual.coefficients
    )
    return coefficients


def pade_candidate(coefficients, numerator_degree, denominator_degree):
    count = len(coefficients)
    if count < numerator_degree+denominator_degree+1:
        return None
    if denominator_degree == 0:
        if any(coefficients[numerator_degree+1:]):
            return None
        return tuple(coefficients[:numerator_degree+1]), (QQ.one(),)
    rows = tuple(range(numerator_degree+1, count))
    matrix_values = Matrix(
        QQ,
        [
            [coefficients[index-j] if index >= j else QQ.zero()
             for j in range(1, denominator_degree+1)]
            for index in rows
        ],
    )
    target = vector(QQ, [-coefficients[index] for index in rows])
    if matrix_values.rank() != denominator_degree:
        return None
    if target not in matrix_values.column_space():
        return None
    denominator = (QQ.one(),)+tuple(matrix_values.solve_right(target))
    numerator = tuple(
        sum(
            denominator[j]*coefficients[index-j]
            for j in range(min(denominator_degree, index)+1)
        )
        for index in range(numerator_degree+1)
    )
    return numerator, denominator


def minimal_pade(coefficients):
    for total_degree in range(2*arguments.max_degree+1):
        candidates = []
        for denominator_degree in range(arguments.max_degree+1):
            numerator_degree = total_degree-denominator_degree
            if 0 <= numerator_degree <= arguments.max_degree:
                result = pade_candidate(
                    coefficients, numerator_degree, denominator_degree
                )
                if result is not None:
                    candidates.append((numerator_degree, denominator_degree, result))
        if len(candidates) == 1:
            return candidates[0]
    return None


h_ring = PolynomialRing(QQ, "h")
h = h_ring.gen()
surface_relation_ring = PolynomialRing(QQ, names=("D", "P", "Q", "E"))
pair_relation_ring = PolynomialRing(QQ, names=("X", "Y"))
pair_X, pair_Y = pair_relation_ring.gens()


def truncated_product(left, right, length):
    result = [QQ.zero() for _ in range(length)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            if i+j < length:
                result[i+j] += left_value*right_value
    return tuple(result)


def surface_relations(surface_series, degree):
    monomials = tuple(
        monomial
        for current_degree in range(degree+1)
        for monomial in surface_relation_ring.monomials_of_degree(current_degree)
    )
    columns = []
    for monomial in monomials:
        series = (QQ.one(),)+(QQ.zero(),)*(arguments.order-1)
        for coordinate, exponent in zip(surface_series, monomial.exponents()[0]):
            for _ in range(exponent):
                series = truncated_product(series, coordinate, arguments.order)
        columns.append(vector(QQ, series))
    relation_kernel = matrix(QQ, columns).transpose().right_kernel_matrix()
    integral_rows = []
    for row in relation_kernel.rows():
        denominator = LCM([value.denominator() for value in row])
        integral = vector(ZZ, [ZZ(denominator*value) for value in row])
        divisor = GCD([abs(value) for value in integral if value])
        integral_rows.append(integral/divisor)
    reduced_rows = matrix(ZZ, integral_rows).LLL().rows() if integral_rows else ()
    relations = tuple(
        primitive_polynomial(sum(
            coefficient*monomial for coefficient, monomial in zip(row, monomials)
        ))
        for row in reduced_rows if row
    )
    return relations


def primitive_polynomial(relation):
    denominators = [coefficient.denominator() for coefficient in relation.coefficients()]
    integral = relation*LCM(denominators)
    coefficient_gcd = GCD([abs(ZZ(value)) for value in integral.coefficients()])
    primitive = integral.parent()(integral/coefficient_gcd)
    if primitive.coefficients()[0] < 0:
        primitive = -primitive
    return primitive


def verify_slope_8_87_mod7_cubic_span(relations):
    """Match the exact centered cubics to the long GF(7) branch certificate."""
    if len(relations) != 2:
        return False
    finite_ring = PolynomialRing(GF(7), names=("D", "P", "Q", "E"))
    D, P, Q, E = finite_ring.gens()
    raw_d, raw_p, raw_q, raw_e = D+3, P+4, Q+3, E+2
    expected = (
        2*raw_d**3 + raw_d**2*raw_p - raw_d*raw_p**2
        - 2*raw_p**3 - raw_d**2*raw_q + 3*raw_d*raw_p*raw_q
        - 3*raw_p**2*raw_q - 3*raw_d*raw_q**2
        - 2*raw_p*raw_q**2 + 3*raw_q**3 - 2*raw_d**2*raw_e
        + raw_d*raw_p*raw_e + 3*raw_d**2 + 3*raw_d*raw_p
        + 3*raw_p**2 - raw_d*raw_q + 2*raw_p*raw_q + raw_q**2
        - raw_d*raw_e - 2*raw_p*raw_e + 2*raw_q*raw_e
        - 2*raw_e**2 - raw_d + 2*raw_p - 3*raw_q + 1,
        3*raw_d**3 - 3*raw_d**2*raw_p - 3*raw_d*raw_p**2
        - raw_p**3 + 2*raw_d**2*raw_q + raw_d*raw_p*raw_q
        + 2*raw_p**2*raw_q - raw_p*raw_q**2 - 2*raw_q**3
        + 3*raw_d**2*raw_e - 3*raw_d*raw_p*raw_e + raw_d**2
        + 3*raw_p**2 - 2*raw_d*raw_q + 3*raw_p*raw_q + raw_q**2
        + 2*raw_d*raw_e - 3*raw_p*raw_e - 2*raw_q*raw_e
        - raw_e**2 - 3*raw_d - 2*raw_p - raw_q + raw_e,
    )
    reduced = []
    for relation in relations:
        image = finite_ring.zero()
        for exponents, coefficient in relation.dict().items():
            image += GF(7)(ZZ(coefficient))*D**exponents[0]*P**exponents[1] \
                *Q**exponents[2]*E**exponents[3]
        reduced.append(image)
    monomials = tuple(
        monomial
        for degree in range(4)
        for monomial in finite_ring.monomials_of_degree(degree)
    )
    rows = lambda polynomials: Matrix(
        GF(7),
        [[polynomial.monomial_coefficient(monomial) for monomial in monomials]
         for polynomial in polynomials],
    )
    expected_rows = rows(expected)
    reduced_rows = rows(reduced)
    return (
        expected_rows.rank() == 2
        and reduced_rows.rank() == 2
        and expected_rows.stack(reduced_rows).rank() == 2
    )


def normalized_relation(row, monomials):
    relation = sum(coefficient*monomial for coefficient, monomial in zip(row, monomials))
    return primitive_polynomial(relation)


def pairwise_relations(surface_series):
    records = []
    for left_index in range(len(surface_series)):
        for right_index in range(left_index+1, len(surface_series)):
            if arguments.pair_only and (
                surface_names[left_index], surface_names[right_index]
            ) != tuple(arguments.pair_only):
                continue
            left = surface_series[left_index]
            right = surface_series[right_index]
            if arguments.pair_total_degree:
                total_degree = arguments.pair_total_degree
                monomials = tuple(
                    pair_X**left_exponent*pair_Y**right_exponent
                    for left_exponent in range(total_degree+1)
                    for right_exponent in range(total_degree+1-left_exponent)
                )
                if len(monomials)+arguments.pair_validation <= arguments.order:
                    columns = []
                    for monomial in monomials:
                        exponents = monomial.exponents()[0]
                        series = (QQ.one(),)+(QQ.zero(),)*(arguments.order-1)
                        for _ in range(exponents[0]):
                            series = truncated_product(series, left, arguments.order)
                        for _ in range(exponents[1]):
                            series = truncated_product(series, right, arguments.order)
                        columns.append(vector(QQ, series))
                    kernel = matrix(QQ, columns).transpose().right_kernel_matrix()
                    for row in kernel.rows():
                        records.append((
                            left_index,
                            right_index,
                            total_degree,
                            total_degree,
                            arguments.order-len(monomials),
                            normalized_relation(row, monomials),
                        ))
                continue
            for left_degree in range(1, arguments.pair_max_degree+1):
                for right_degree in range(1, arguments.pair_max_degree+1):
                    monomials = tuple(
                        pair_X**left_exponent*pair_Y**right_exponent
                        for left_exponent in range(left_degree+1)
                        for right_exponent in range(right_degree+1)
                    )
                    if len(monomials)+arguments.pair_validation > arguments.order:
                        continue
                    columns = []
                    for monomial in monomials:
                        exponents = monomial.exponents()[0]
                        series = (QQ.one(),)+(QQ.zero(),)*(arguments.order-1)
                        for _ in range(exponents[0]):
                            series = truncated_product(series, left, arguments.order)
                        for _ in range(exponents[1]):
                            series = truncated_product(series, right, arguments.order)
                        columns.append(vector(QQ, series))
                    kernel = matrix(QQ, columns).transpose().right_kernel_matrix()
                    for row in kernel.rows():
                        records.append((
                            left_index,
                            right_index,
                            left_degree,
                            right_degree,
                            arguments.order-len(monomials),
                            normalized_relation(row, monomials),
                        ))
    return tuple(records)


def supported_plane_relation(surface_series):
    """Recover an exact relation on the support of a certified GF(7) eliminant."""
    artifact = json.loads(Path(arguments.supported_plane_ideal).read_text())
    if artifact.get("schema") != "q80-cm24-formal-branch-ideal-v1":
        raise ValueError("unexpected supported-plane ideal schema")
    finite_ring = PolynomialRing(
        GF(7), names=tuple(artifact["affine_variables"]), order="degrevlex"
    )
    finite_ideal = finite_ring.ideal(tuple(
        finite_ring(polynomial) for polynomial in artifact["affine_generators"]
    ))
    keep_names = tuple(arguments.supported_plane_pair)
    if len(set(keep_names)) != 2 or not set(keep_names) <= set(surface_names):
        raise ValueError("supported plane pair must contain two distinct surface names")
    eliminated = tuple(
        variable for variable in finite_ring.gens() if str(variable) not in keep_names
    )
    elimination_ideal = finite_ideal.elimination_ideal(eliminated)
    if len(elimination_ideal.gens()) != 1:
        raise ValueError("supported plane elimination is not principal")
    modular_polynomial = elimination_ideal.gens()[0]
    support = tuple(sorted(modular_polynomial.dict()))
    if arguments.order <= len(support):
        raise ValueError(
            f"supported relation needs more than {len(support)} exact coefficients"
        )

    exact_pair_ring = PolynomialRing(QQ, names=keep_names)
    exact_variables = exact_pair_ring.gens()
    exact_monomials = tuple(
        exact_variables[0]**exponents[finite_ring.variable_names().index(keep_names[0])]
        * exact_variables[1]**exponents[finite_ring.variable_names().index(keep_names[1])]
        for exponents in support
    )
    pair_series = tuple(surface_series[surface_names.index(name)] for name in keep_names)
    columns = []
    for monomial in exact_monomials:
        series = (QQ.one(),)+(QQ.zero(),)*(arguments.order-1)
        for coordinate, exponent in zip(pair_series, monomial.exponents()[0]):
            for _ in range(exponent):
                series = truncated_product(series, coordinate, arguments.order)
        columns.append(vector(QQ, series))
    kernel = Matrix(QQ, columns).transpose().right_kernel_matrix()
    if kernel.nrows() != 1:
        return None, len(support), arguments.order-len(support), kernel.nrows()
    exact_relation = normalized_relation(kernel.row(0), exact_monomials)

    modular_pair_ring = PolynomialRing(GF(7), names=keep_names)
    modular_image = modular_pair_ring.zero()
    for exponents, coefficient in exact_relation.dict().items():
        modular_image += (
            GF(7)(ZZ(coefficient))
            * modular_pair_ring.gen(0)**exponents[0]
            * modular_pair_ring.gen(1)**exponents[1]
        )
    expected_image = modular_pair_ring(str(modular_polynomial))
    quotient = next(
        modular_image.monomial_coefficient(monomial)
        / expected_image.monomial_coefficient(monomial)
        for monomial in expected_image.monomials()
        if expected_image.monomial_coefficient(monomial)
    )
    if modular_image != quotient*expected_image:
        raise AssertionError("exact supported relation does not reduce to GF(7) eliminant")
    return exact_relation, len(support), arguments.order-len(support), kernel.nrows()


selected_slopes = tuple(map(QQ, arguments.slope)) if arguments.slope else (QQ(8)/87, QQ(1)/12)
surface_series_records = []
for slope in selected_slopes:
    coefficients = lift_branch_fast(slope)
    surface_series = []
    for column in surface_columns:
        values = tuple(coefficient[column] for coefficient in coefficients)
        rational_values = tuple(value.rational() for value in values)
        surface_series.append(rational_values)
    surface_centers = (QQ(-1)/2, QQ(9)/4, QQ(-9)/4, QQ(-27)/32)
    centered_surface_series = []
    for values, center in zip(surface_series, surface_centers):
        centered_surface_series.append((values[0]-center,)+values[1:])
    surface_series_records.append(
        {
            "slope": str(slope),
            "order": arguments.order,
            "normalization": "P=9/4+h",
            "centers": {
                name: str(center)
                for name, center in zip(surface_names, surface_centers)
            },
            "series": {
                name: [str(value) for value in values]
                for name, values in zip(surface_names, surface_series)
            },
        }
    )
    if not arguments.quiet_series:
        print(
            f"Q80CM24BRANCH|slope={slope}|order={arguments.order}|"
            +"|".join(
                f"{name}_series={','.join(map(str, values))}"
                for name, values in zip(surface_names, surface_series)
            ),
            flush=True,
        )
    if arguments.surface_series_only:
        print(
            f"Q80CM24BRANCH|slope={slope}|order={arguments.order}|"
            "surface_series_export_ready=1|diagnostics=SKIPPED",
            flush=True,
        )
        continue
    for relation_degree in range(1, arguments.relation_degree+1):
        relations = surface_relations(centered_surface_series, relation_degree)
        monomial_count = binomial(4+relation_degree, relation_degree)
        determined = arguments.order > monomial_count
        print(
            f"Q80CM24BRANCH|slope={slope}|relation_degree={relation_degree}|"
            f"monomials={monomial_count}|determined={int(determined)}|"
            f"coordinates=centered_at_CM24|kernel={len(relations)}|"
            f"relations={relations if determined else ()}",
            flush=True,
        )
        if slope == QQ(8)/87 and relation_degree == 3 and determined:
            mod7_span_match = verify_slope_8_87_mod7_cubic_span(relations)
            assert mod7_span_match
            print(
                "Q80CM24BRANCH|slope=8/87|relation_degree=3|"
                "mod7_order145_reference_span_match=1",
                flush=True,
            )
    if arguments.supported_plane_ideal:
        relation, support_size, validation, kernel_dimension = supported_plane_relation(
            centered_surface_series
        )
        print(
            f"Q80CM24PLANE|slope={slope}|pair="
            f"{','.join(arguments.supported_plane_pair)}|support={support_size}|"
            f"validation={validation}|kernel={kernel_dimension}|relation={relation}|"
            f"status={'PASS' if relation is not None else 'NO_RELATION_ON_MOD7_SUPPORT'}",
            flush=True,
        )
    pair_records = pairwise_relations(centered_surface_series)
    print(
        f"Q80CM24BRANCH|slope={slope}|pair_max_degree={arguments.pair_max_degree}|"
        f"pair_validation={arguments.pair_validation}|pair_candidates={len(pair_records)}",
        flush=True,
    )
    for left_index, right_index, left_degree, right_degree, validation, relation in pair_records:
        print(
            f"Q80CM24BRANCH|slope={slope}|pair="
            f"{surface_names[left_index]},{surface_names[right_index]}|"
            f"bidegree={left_degree},{right_degree}|validation={validation}|"
            f"relation={relation}",
            flush=True,
        )
    for name, values in zip(surface_names, surface_series):
        candidate = minimal_pade(values)
        if candidate is None:
            continue
        numerator_degree, denominator_degree, (numerator, denominator) = candidate
        numerator_polynomial = sum(value*h**index for index, value in enumerate(numerator))
        denominator_polynomial = sum(value*h**index for index, value in enumerate(denominator))
        print(
            f"Q80CM24BRANCH|slope={slope}|coordinate={name}|"
            f"pade_degrees={numerator_degree},{denominator_degree}|"
            f"pade=({numerator_polynomial})/({denominator_polynomial})",
            flush=True,
        )

if arguments.surface_series_output:
    output_path = Path(arguments.surface_series_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(
            {
                "schema": "q80-cm24-qq-surface-series-v1",
                "scope": "exact_formal_series_not_global_algebraization",
                "records": surface_series_records,
            },
            indent=2,
            sort_keys=True,
        )+"\n"
    )
    print(
        f"Q80CM24BRANCH|surface_series_output={output_path}|"
        f"records={len(surface_series_records)}",
        flush=True,
    )

Jet2.precision = 3
print(
    f"Q80CM24BRANCH|branches={len(selected_slopes)}|order={arguments.order}|"
    "normalization=P-9/4|status=PASS",
    flush=True,
)
