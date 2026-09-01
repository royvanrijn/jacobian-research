#!/usr/bin/env sage -python
"""Compile one p-adic Riemann--Roch sample in the exact quadratic gauge.

The exact base value U=16+7*omega reduces, under the pinned p=19 legacy
PGL2 alignment, to the already certified local sample T=1.  The worker uses
the generic integral basis (1,z,e) modulo 19^5, lifts the simple and double
infinity branches, and constructs normalized generators of L(2P),L(3P).
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    EllipticCurve,
    GF,
    Matrix,
    PolynomialRing,
    Qp,
    FunctionField,
    LaurentSeriesRing,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_SOURCE = RESULTS / "q80-third-q12-exact-pencil-p19-adic-precision64.json"
DEFAULT_LIFT = RESULTS / "q80-third-q12-discriminant-factors-p19-adic-precision5.json"
DEFAULT_BASIS = RESULTS / "q80-third-q12-integral-basis-mod19-power.json"
DEFAULT_CONTROL = RESULTS / "q80-third-q12-um2-p19-weierstrass-sample.json"
DEFAULT_OUTPUT = RESULTS / "q80-third-q12-riemann-roch-p19-adic-sample.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--lift", type=Path, default=DEFAULT_LIFT)
parser.add_argument("--basis", type=Path, default=DEFAULT_BASIS)
parser.add_argument("--control", type=Path, default=DEFAULT_CONTROL)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--base-constant", type=int, default=16)
parser.add_argument("--base-anti", type=int, default=7)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
for name in ("source", "lift", "basis", "control", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


source = json.loads(args.source.read_text())
lift = json.loads(args.lift.read_text())
basis_certificate = json.loads(args.basis.read_text())
control = json.loads(args.control.read_text())
if source.get("status") != "PASS_EXACT_THIRD_Q12_PENCIL_REDUCTION_MOD_19_POWER":
    raise ValueError("exact p-adic pencil is not certified")
if lift.get("status") != "PASS_EXACT_THIRD_Q12_DISCRIMINANT_FACTOR_HENSEL_LIFT_P19":
    raise ValueError("conductor lift is not certified")
if basis_certificate.get("status") != "PASS_EXACT_THIRD_Q12_GENERIC_INTEGRAL_BASIS_MOD19_POWER":
    raise ValueError("generic integral basis is not certified")
if control.get("status") != "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_MOD19_QUADRATIC":
    raise ValueError("legacy finite-field Weierstrass sample is not certified")
if control["specialization"]["new_base_coefficients_1_r"] != [1, 0]:
    raise ValueError("legacy positive control is not the T=1 sample")
if basis_certificate["inputs"]["source"]["sha256"] != sha256(args.source):
    raise ValueError("integral-basis source hash changed")
if basis_certificate["inputs"]["lift"]["sha256"] != sha256(args.lift):
    raise ValueError("integral-basis lift hash changed")

prime = 19
digits = int(lift["specialization"]["digits"])
modulus = prime**digits
padic = Qp(prime, prec=digits)
omega_ring = PolynomialRing(padic, "omega_polynomial")
omega_polynomial = omega_ring.gen()
omega_square = int(source["quadratic_field"]["omega_square_modulus"])
quadratic = padic.extension(
    omega_polynomial**2 - padic(omega_square), names="omega"
)
omega = quadratic.gen()


def quadratic_element(coordinates):
    return quadratic(int(coordinates[0])) + quadratic(int(coordinates[1])) * omega


def residue_coordinates(value):
    coefficients = list(quadratic(value).polynomial().list())
    coefficients += [padic.zero()] * (2 - len(coefficients))
    precisions = []
    for coefficient in coefficients:
        try:
            precisions.append(min(digits, int(coefficient.precision_absolute())))
        except TypeError:
            precisions.append(digits)
    precision = min(precisions)
    if precision < digits:
        raise ArithmeticError(f"lost p-adic precision: only {precision} digits remain")
    return [int(coefficient.lift()) % modulus for coefficient in coefficients]


# The default is the aligned legacy T=1 value.  Other residue classes exercise
# the same generic p-adic producer and derive their infinity branches below.
base_constant = int(args.base_constant) % modulus
base_anti = int(args.base_anti) % modulus
U_value = quadratic(base_constant) + quadratic(base_anti) * omega


def evaluate_u_polynomial(record):
    return sum(
        quadratic_element(coefficient) * U_value**index
        for index, coefficient in enumerate(record)
    )


def evaluate_rational_record(record):
    numerator = evaluate_u_polynomial(
        record["numerator_coefficients_low_to_high_U_1_omega"]
    )
    denominator = evaluate_u_polynomial(
        record["denominator_coefficients_low_to_high_U_1_omega"]
    )
    if denominator.valuation() != 0:
        raise ZeroDivisionError("specialized U denominator is not a unit")
    return numerator / denominator


old_function = FunctionField(quadratic, "W")
W = old_function.gen()
x_ring = PolynomialRing(old_function, "x")
x = x_ring.gen()
equation = x_ring.zero()
for u_degree, w_degree, x_degree, coordinates in source["pencil"][
    "terms_V_W_old_x_coefficient_1_omega"
]:
    equation += (
        quadratic_element(coordinates)
        * U_value**u_degree
        * W**w_degree
        * x**x_degree
    )
if equation.degree() != 3:
    raise ArithmeticError("specialized exact pencil is not cubic")
monic_equation = equation.monic()


def specialized_w_polynomial(record):
    return sum(
        old_function(evaluate_rational_record(coefficient)) * W**index
        for index, coefficient in enumerate(record["coefficients_low_to_high_W"])
    )


L = specialized_w_polynomial(lift["factorization"]["L"])
Q = specialized_w_polynomial(lift["factorization"]["Q"])
A = specialized_w_polynomial(lift["integral_basis_candidate"]["A"])
B = specialized_w_polynomial(lift["integral_basis_candidate"]["B"])
conductor = L * Q

# Only vector-space operations are needed for L(0P),...,L(3P).  Avoid Sage's
# function-field constructor, which tries an unavailable p-adic factorization,
# by storing elements directly in the old basis (1,z,z^2).
ZERO_E = (old_function.zero(),) * 3
ONE_E = (old_function.one(), old_function.zero(), old_function.zero())
Z_E = (old_function.zero(), old_function.one(), old_function.zero())


def e_add(left, right):
    return tuple(left[index] + right[index] for index in range(3))


def e_scale(value, scalar):
    return tuple(scalar * coefficient for coefficient in value)


e = (B / conductor, A / conductor, old_function.one() / conductor)
integral_basis = [ONE_E, Z_E, e]

# Weighted infinity chart t=1/W, xi=t^3*z.
infinity_ring = PolynomialRing(quadratic, names=("t", "xi"))
t, xi = infinity_ring.gens()
infinity_equation = infinity_ring.zero()
for u_degree, w_degree, x_degree, coordinates in source["pencil"][
    "terms_V_W_old_x_coefficient_1_omega"
]:
    infinity_equation += (
        quadratic_element(coordinates)
        * U_value**u_degree
        * t ** (9 - w_degree - 3 * x_degree)
        * xi**x_degree
    )

# Derive the distinguished simple branch and the conjugate pair above the
# double branch from the weighted infinity cubic modulo 19.
finite_base = GF(prime)
finite_modulus_ring = PolynomialRing(finite_base, "omega_residue_polynomial")
omega_residue_polynomial = finite_modulus_ring.gen()
finite = GF(
    prime**2,
    "omega_residue",
    modulus=omega_residue_polynomial**2 - finite_base(omega_square % prime),
)
omega_residue = finite.gen()


def finite_reduction(value):
    constant, anti = residue_coordinates(value)
    return finite(constant % prime) + finite(anti % prime) * omega_residue


infinity_residue_ring = PolynomialRing(finite, "xi_residue")
xi_residue = infinity_residue_ring.gen()
infinity_cubic_residue = sum(
    finite_reduction(coefficient) * xi_residue**exponent[1]
    for exponent, coefficient in infinity_equation(t=0).dict().items()
)
infinity_roots = infinity_cubic_residue.roots()
simple_roots = [root for root, multiplicity in infinity_roots if multiplicity == 1]
double_roots = [root for root, multiplicity in infinity_roots if multiplicity == 2]
if len(simple_roots) != 1 or len(double_roots) != 1:
    raise ArithmeticError(
        f"weighted infinity cubic does not have one simple and one double root: {infinity_roots}"
    )


def lift_finite(value):
    coefficients = list(finite(value).list())
    coefficients += [finite_base.zero()] * (2 - len(coefficients))
    return quadratic(int(coefficients[0])) + quadratic(int(coefficients[1])) * omega


simple_seed = lift_finite(simple_roots[0])
double_seed = lift_finite(double_roots[0])

series_precision = 45
simple_series_ring = LaurentSeriesRing(quadratic, "s", default_prec=series_precision)
s = simple_series_ring.gen()
simple_x_ring = PolynomialRing(simple_series_ring, "X")
X = simple_x_ring.gen()
simple_equation = simple_x_ring(infinity_equation(s, X))
simple_root = simple_series_ring(simple_seed)
for unused in range(8):
    simple_root -= simple_equation(simple_root) / simple_equation.derivative()(simple_root)
if simple_equation(simple_root).valuation() < series_precision - 6:
    raise ArithmeticError("simple infinity branch did not converge")

# At the double branch use t=s^2 and xi=xi0+c2*s^2+c3*s^3+....
branch_parameter_ring = PolynomialRing(quadratic, names=("c2", "c3"))
c2_variable, c3_variable = branch_parameter_ring.gens()
trial_ring = PolynomialRing(branch_parameter_ring, "s_trial")
s_trial = trial_ring.gen()
trial = trial_ring(
    infinity_equation(
        s_trial**2,
        double_seed + c2_variable * s_trial**2 + c3_variable * s_trial**3,
    )
)
c2_ring = PolynomialRing(quadratic, "c2")
c2_polynomial = c2_ring(trial[4](c2_variable, 0))
c2_roots = c2_polynomial.roots(multiplicities=False)
if len(c2_roots) != 1:
    raise ArithmeticError("double branch has no unique p-adic quadratic shift")
c2_value = c2_roots[0]
c3_ring = PolynomialRing(quadratic, "c3")
c3_polynomial = None
for exponent in range(5, 10):
    candidate = c3_ring(trial[exponent](c2_value, c3_variable))
    if candidate:
        c3_polynomial = candidate.monic()
        break
if c3_polynomial is None or c3_polynomial.degree() != 2:
    raise ArithmeticError("double branch has no quadratic p-adic cubic shift")
c3_roots = c3_polynomial.roots(multiplicities=False)
if c3_roots:
    branch_field = quadratic
    c3_value = c3_roots[0]
    branch_degree = 1
    branch_map = quadratic
    c3_values = [c3_value]
else:
    # Sage does not implement this nested q-adic extension directly.  Flatten
    # it to the norm quartic over Q_19, then recover the embedded omega from
    # the quadratic c3 equation.
    def quadratic_parts(value):
        coefficients = list(quadratic(value).polynomial().list())
        coefficients += [padic.zero()] * (2 - len(coefficients))
        return coefficients[0], coefficients[1]

    def quadratic_conjugate(value):
        constant, anti = quadratic_parts(value)
        return quadratic(constant) - quadratic(anti) * omega

    conjugate_c3_polynomial = sum(
        quadratic_conjugate(coefficient) * c3_ring.gen()**index
        for index, coefficient in enumerate(c3_polynomial.list())
    )
    norm_polynomial_quadratic = c3_polynomial * conjugate_c3_polynomial
    norm_ring = PolynomialRing(padic, "c3_norm")
    norm_polynomial = norm_ring.zero()
    for index, coefficient in enumerate(norm_polynomial_quadratic.list()):
        constant, anti = quadratic_parts(coefficient)
        if anti.valuation() < digits:
            raise ArithmeticError("c3 norm polynomial did not descend to Q_19")
        norm_polynomial += constant * norm_ring.gen()**index
    norm_polynomial = norm_polynomial.monic()
    branch_field = padic.extension(norm_polynomial, names="c3_branch")
    c3_value = branch_field.gen()
    constant_part = branch_field.zero()
    anti_part = branch_field.zero()
    for index, coefficient in enumerate(c3_polynomial.list()):
        constant, anti = quadratic_parts(coefficient)
        constant_part += branch_field(constant) * c3_value**index
        anti_part += branch_field(anti) * c3_value**index
    if not anti_part or anti_part.valuation() != 0:
        raise ArithmeticError("cannot recover omega from flattened c3 extension")
    omega_image = -constant_part / anti_part
    if (omega_image**2 - branch_field(omega_square)).valuation() < digits:
        raise ArithmeticError("flattened c3 extension has the wrong omega embedding")

    def branch_map(value):
        constant, anti = quadratic_parts(value)
        return branch_field(constant) + branch_field(anti) * omega_image

    c3_values = [
        c3_value,
        -branch_map(c3_polynomial[1]) - c3_value,
    ]
    branch_degree = 2

double_polynomial_ring = PolynomialRing(branch_field, names=("td", "xid"))
td, xid = double_polynomial_ring.gens()
double_equation_polynomial = double_polynomial_ring.zero()
for exponent, coefficient in infinity_equation.dict().items():
    double_equation_polynomial += (
        branch_map(coefficient) * td**exponent[0] * xid**exponent[1]
    )
double_branch_data = []
for branch_index, c3_branch_value in enumerate(c3_values):
    double_series_ring = LaurentSeriesRing(
        branch_field, f"sd{branch_index}", default_prec=2 * series_precision
    )
    sd = double_series_ring.gen()
    double_x_ring = PolynomialRing(double_series_ring, f"Xd{branch_index}")
    Xd = double_x_ring.gen()
    double_equation = double_x_ring(double_equation_polynomial(sd**2, Xd))
    double_root = double_series_ring(
        branch_map(double_seed) + branch_map(c2_value) * sd**2 + c3_branch_value * sd**3
    )
    for unused in range(9):
        double_root -= double_equation(double_root) / double_equation.derivative()(double_root)
    if double_equation(double_root).valuation() < 2 * series_precision - 9:
        raise ArithmeticError("double infinity branch did not converge")
    double_branch_data.append((sd, double_root))


def evaluate_polynomial(polynomial, w_series, mapper):
    return sum(
        w_series.parent()(mapper(coefficient)) * w_series**index
        for index, coefficient in enumerate(polynomial.list())
    )


def evaluate_rational(value, w_series, mapper):
    return evaluate_polynomial(value.numerator(), w_series, mapper) / evaluate_polynomial(
        value.denominator(), w_series, mapper
    )


def evaluate_element(value, w_series, z_series, mapper):
    return sum(
        evaluate_rational(coefficient, w_series, mapper) * z_series**index
        for index, coefficient in enumerate(value)
    )


simple_basis = [
    evaluate_element(value, s**-1, simple_root * s**-3, lambda coefficient: coefficient)
    for value in integral_basis
]
double_bases = [
    [
        evaluate_element(value, sd**-2, double_root * sd**-6, branch_map)
        for value in integral_basis
    ]
    for sd, double_root in double_branch_data
]
integral_basis_valuations = [
    [int(simple_basis[index].valuation()), int(double_bases[0][index].valuation())]
    for index in range(3)
]
if integral_basis_valuations != [[0, 0], [-3, -6], [-2, -4]]:
    raise ArithmeticError(
        f"unexpected infinity valuations of p-adic integral basis: {integral_basis_valuations}"
    )

degree_bound = 3
labels = [
    (basis_index, degree)
    for basis_index in range(3)
    for degree in range(degree_bound + 1)
]
functions = [
    e_scale(integral_basis[basis_index], W**degree)
    for basis_index, degree in labels
]
simple_columns = [
    simple_basis[basis_index] * s**-degree
    for basis_index, degree in labels
]
double_column_sets = [
    [
        double_basis[basis_index] * sd ** (-2 * degree)
        for basis_index, degree in labels
    ]
    for (sd, unused), double_basis in zip(double_branch_data, double_bases)
]


def riemann_roch_space(pole_bound):
    rows = []
    for exponent in range(
        min(value.valuation() for value in simple_columns), -pole_bound
    ):
        rows.append([branch_map(value[exponent]) for value in simple_columns])
    for double_columns in double_column_sets:
        for exponent in range(min(value.valuation() for value in double_columns), 0):
            rows.append([value[exponent] for value in double_columns])
    return Matrix(branch_field, rows).right_kernel()


spaces = [riemann_roch_space(bound) for bound in range(4)]
dimensions = [int(space.dimension()) for space in spaces]
if dimensions != [1, 1, 2, 3]:
    raise ArithmeticError(f"unexpected p-adic Riemann--Roch dimensions {dimensions}")


simple_branch_series_ring = LaurentSeriesRing(
    branch_field, "simple_branch_s", default_prec=series_precision
)
simple_branch_s = simple_branch_series_ring.gen()


def map_simple_series(value):
    if not value:
        return simple_branch_series_ring.zero()
    try:
        upper = min(series_precision, int(value.prec()))
    except TypeError:
        upper = int(value.degree()) + 1
    return sum(
        (
            simple_branch_series_ring(branch_map(value[exponent]))
            * simple_branch_s**exponent
            for exponent in range(int(value.valuation()), upper)
        ),
        simple_branch_series_ring.zero(),
    )


simple_columns_over_branch = [map_simple_series(value) for value in simple_columns]


def simple_expansion(coordinates):
    return sum(
        (
            simple_branch_series_ring(coefficient) * column
            for coefficient, column in zip(coordinates, simple_columns_over_branch)
        ),
        simple_branch_series_ring.zero(),
    )


def simple_valuation(value):
    return int(simple_expansion(value).valuation())


x_coordinates = next(
    coordinates
    for coordinates in spaces[2].basis()
    if simple_valuation(coordinates) == -2
)
y_coordinates = next(
    coordinates
    for coordinates in spaces[3].basis()
    if simple_valuation(coordinates) == -3
)
x_coordinates = vector(branch_field, x_coordinates)
y_coordinates = vector(branch_field, y_coordinates)
one_coordinates = vector(
    branch_field,
    [1 if label == (0, 0) else 0 for label in labels],
)

x_local = simple_expansion(x_coordinates)
x_coordinates /= x_local[-2]
x_local = simple_expansion(x_coordinates)
x_coordinates -= x_local[0] * one_coordinates
x_local = simple_expansion(x_coordinates)
if x_local[-2] != 1 or x_local[0].valuation() < digits:
    raise ArithmeticError("failed to normalize p-adic pole-two generator")

y_local = simple_expansion(y_coordinates)
y_coordinates /= y_local[-3]
y_local = simple_expansion(y_coordinates)
y_coordinates -= y_local[-2] * x_coordinates
y_local = simple_expansion(y_coordinates)
y_coordinates -= y_local[0] * one_coordinates
y_local = simple_expansion(y_coordinates)
if y_local[-3] != 1 or y_local[-2].valuation() < digits or y_local[0].valuation() < digits:
    raise ArithmeticError("failed to normalize p-adic pole-three generator")


def quadratic_preimage(value):
    if branch_degree == 1:
        return quadratic(value)
    value_coefficients = list(branch_field(value).polynomial().list())
    omega_coefficients = list(branch_field(omega_image).polynomial().list())
    degree = int(branch_field.degree())
    value_coefficients += [padic.zero()] * (degree - len(value_coefficients))
    omega_coefficients += [padic.zero()] * (degree - len(omega_coefficients))
    pivot = next(
        index
        for index in range(1, degree)
        if omega_coefficients[index].valuation() == 0
    )
    anti = value_coefficients[pivot] / omega_coefficients[pivot]
    constant = value_coefficients[0] - anti * omega_coefficients[0]
    residual = branch_field(value) - branch_field(constant) - branch_field(anti) * omega_image
    if residual.valuation() < digits:
        raise ArithmeticError(
            "normalized Riemann--Roch coordinate did not descend: "
            f"residual valuation {residual.valuation()}"
        )
    return quadratic(constant) + quadratic(anti) * omega


x_quadratic_coordinates = [quadratic_preimage(value) for value in x_coordinates]
y_quadratic_coordinates = [quadratic_preimage(value) for value in y_coordinates]


def element_from_coordinates(coordinates):
    result = ZERO_E
    for coefficient, function in zip(coordinates, functions):
        result = e_add(result, e_scale(function, coefficient))
    return result


def e_mul(left, right):
    coefficients = [old_function.zero()] * 5
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            coefficients[left_degree + right_degree] += left_coefficient * right_coefficient
    # z^3=-b*z^2-c*z-d for the monic cubic z^3+b*z^2+c*z+d.
    b_cubic = monic_equation[2]
    c_cubic = monic_equation[1]
    d_cubic = monic_equation[0]
    for degree in (4, 3):
        coefficient = coefficients[degree]
        coefficients[degree] = old_function.zero()
        coefficients[degree - 3] -= coefficient * d_cubic
        coefficients[degree - 2] -= coefficient * c_cubic
        coefficients[degree - 1] -= coefficient * b_cubic
    return tuple(coefficients[:3])


def e_pow(value, exponent):
    result = ONE_E
    base = value
    while exponent:
        if exponent & 1:
            result = e_mul(result, base)
        base = e_mul(base, base)
        exponent >>= 1
    return result


x_element = element_from_coordinates(x_quadratic_coordinates)
y_element = element_from_coordinates(y_quadratic_coordinates)
relation_terms = (
    e_pow(y_element, 2),
    e_mul(x_element, y_element),
    y_element,
    e_scale(e_pow(x_element, 3), -1),
    e_scale(e_pow(x_element, 2), -1),
    e_scale(x_element, -1),
    e_scale(ONE_E, -1),
)


def constant_relation_kernel(function_values):
    rows = []
    for component in range(3):
        values = [value[component] for value in function_values]
        numerators = [value.numerator() for value in values]
        denominators = [value.denominator() for value in values]
        cleared_numerators = []
        for index, numerator in enumerate(numerators):
            cleared = numerator
            for other_index, denominator in enumerate(denominators):
                if other_index != index:
                    cleared *= denominator
            cleared_numerators.append(cleared)
        maximum_degree = max(
            (-1 if not numerator else int(numerator.degree()))
            for numerator in cleared_numerators
        )
        for degree in range(maximum_degree + 1):
            rows.append([numerator[degree] for numerator in cleared_numerators])
    return Matrix(quadratic, rows).right_kernel()


relation_kernel = constant_relation_kernel(relation_terms)
if relation_kernel.dimension() != 1:
    raise ArithmeticError(
        f"p-adic Weierstrass relation has dimension {relation_kernel.dimension()}"
    )
relation = vector(quadratic, relation_kernel.basis()[0])
if relation[0].valuation() != 0 or relation[3].valuation() != 0:
    raise ArithmeticError("p-adic Weierstrass relation has non-unit leading terms")
relation /= relation[0]
cubic_scale = relation[3]
a1 = relation[1]
a2 = relation[4]
a3 = relation[2] * cubic_scale
a4 = relation[5] * cubic_scale
a6 = relation[6] * cubic_scale**2
curve = EllipticCurve(quadratic, [a1, a2, a3, a4, a6])
if curve.discriminant().valuation() != 0:
    raise ArithmeticError("p-adic Weierstrass sample is not smooth")

weierstrass_x = e_scale(x_element, cubic_scale)
weierstrass_y = e_scale(y_element, cubic_scale)
literal_relation = e_add(
    e_add(
        e_add(e_pow(weierstrass_y, 2), e_scale(e_mul(weierstrass_x, weierstrass_y), a1)),
        e_scale(weierstrass_y, a3),
    ),
    e_scale(
        e_add(
            e_add(
                e_add(e_pow(weierstrass_x, 3), e_scale(e_pow(weierstrass_x, 2), a2)),
                e_scale(weierstrass_x, a4),
            ),
            e_scale(ONE_E, a6),
        ),
        -1,
    ),
)


def gauss_valuation(polynomial):
    if not polynomial:
        return digits
    return min(coefficient.valuation() for coefficient in polynomial.list())


for coefficient in literal_relation:
    if gauss_valuation(coefficient.denominator()) != 0:
        raise ArithmeticError("literal-replay denominator is not a p-adic unit")
    if gauss_valuation(coefficient.numerator()) < digits:
        raise ArithmeticError("literal Weierstrass replay failed modulo 19^5")


def local_r_coordinates_mod19(value):
    constant, anti = residue_coordinates(value)
    # omega -> 16*(2*r+12)=2+13*r modulo 19.
    return [(constant + 2 * anti) % prime, (13 * anti) % prime]


weierstrass_coefficients = [a1, a2, a3, a4, a6]
local_weierstrass = [local_r_coordinates_mod19(value) for value in weierstrass_coefficients]
is_legacy_sample = base_constant % prime == 16 and base_anti % prime == 7
if is_legacy_sample and local_weierstrass != control["weierstrass"]["a1_a2_a3_a4_a6"]:
    raise ArithmeticError(
        "p-adic Weierstrass equation does not reduce to the legacy T=1 positive control"
    )


def rational_function_record(value):
    value = old_function(value)
    return {
        "numerator_coefficients_low_to_high_W_mod_19_power_1_omega": [
            residue_coordinates(coefficient) for coefficient in value.numerator().list()
        ],
        "denominator_coefficients_low_to_high_W_mod_19_power_1_omega": [
            residue_coordinates(coefficient) for coefficient in value.denominator().list()
        ],
    }


def function_record(value):
    return {
        "coefficients_low_to_high_z": [
            rational_function_record(coefficient) for coefficient in value
        ]
    }


conductor_polynomial = conductor.numerator()
if conductor.denominator() != 1 or conductor_polynomial.degree() != 5:
    raise ArithmeticError("specialized conductor is not a degree-five polynomial")


def polynomial_zero_mod_19_power(value):
    return not value or gauss_valuation(value) >= digits


def rational_element_to_conductor_power(value, exponent):
    expected_denominator = conductor_polynomial**exponent
    result = []
    for coefficient in value:
        if not coefficient:
            result.append(expected_denominator.parent().zero())
            continue
        denominator = coefficient.denominator()
        scalar = denominator.leading_coefficient() / expected_denominator.leading_coefficient()
        difference = denominator - scalar * expected_denominator
        if not polynomial_zero_mod_19_power(difference):
            raise ArithmeticError(
                f"forward-map denominator is not a scalar times (LQ)^{exponent}"
            )
        result.append(coefficient.numerator() / scalar)
    return exponent, tuple(result)


X_H = rational_element_to_conductor_power(weierstrass_x, 2)
Y_H = rational_element_to_conductor_power(weierstrass_y, 3)
ONE_H = (0, (conductor_polynomial.parent().one(),) + (conductor_polynomial.parent().zero(),) * 2)
W_H = (0, (conductor_polynomial.parent().gen(),) + (conductor_polynomial.parent().zero(),) * 2)
Z_H = (0, (conductor_polynomial.parent().zero(), conductor_polynomial.parent().one(), conductor_polynomial.parent().zero()))
cubic_polynomials = []
for coefficient in (monic_equation[2], monic_equation[1], monic_equation[0]):
    if coefficient.denominator() != 1:
        raise ArithmeticError("monic cubic coefficient is not polynomial in W")
    cubic_polynomials.append(coefficient.numerator())
b_cubic_polynomial, c_cubic_polynomial, d_cubic_polynomial = cubic_polynomials


def h_scale(value, scalar):
    return value[0], tuple(scalar * coefficient for coefficient in value[1])


def h_add(left, right):
    exponent = max(left[0], right[0])
    left_factor = conductor_polynomial ** (exponent - left[0])
    right_factor = conductor_polynomial ** (exponent - right[0])
    return exponent, tuple(
        left[1][index] * left_factor + right[1][index] * right_factor
        for index in range(3)
    )


def h_mul(left, right):
    coefficients = [conductor_polynomial.parent().zero()] * 5
    for left_degree, left_coefficient in enumerate(left[1]):
        for right_degree, right_coefficient in enumerate(right[1]):
            coefficients[left_degree + right_degree] += left_coefficient * right_coefficient
    for degree in (4, 3):
        coefficient = coefficients[degree]
        coefficients[degree] = conductor_polynomial.parent().zero()
        coefficients[degree - 3] -= coefficient * d_cubic_polynomial
        coefficients[degree - 2] -= coefficient * c_cubic_polynomial
        coefficients[degree - 1] -= coefficient * b_cubic_polynomial
    return left[0] + right[0], tuple(coefficients[:3])


def h_pow(value, exponent):
    result = ONE_H
    base = value
    while exponent:
        if exponent & 1:
            result = h_mul(result, base)
        base = h_mul(base, base)
        exponent >>= 1
    return result


def h_evaluate(value, w_value):
    denominator = conductor_polynomial(w_value) ** value[0]
    if denominator.valuation() != 0:
        raise ZeroDivisionError
    return [coefficient(w_value) / denominator for coefficient in value[1]]


def h_is_zero_mod_19_power(value):
    return all(polynomial_zero_mod_19_power(coefficient) for coefficient in value[1])


def h_relation_kernel(function_values):
    rows = []
    good_values = []
    for anti in range(prime):
        for constant in range(prime):
            w_value = quadratic(constant) + quadratic(anti) * omega
            if conductor_polynomial(w_value).valuation() != 0:
                continue
            evaluated = [h_evaluate(value, w_value) for value in function_values]
            good_values.append([constant, anti])
            for component in range(3):
                rows.append([value[component] for value in evaluated])
            if len(good_values) >= 24:
                kernel = Matrix(quadratic, rows).right_kernel()
                if kernel.dimension() <= 1:
                    return kernel, good_values
    return Matrix(quadratic, rows).right_kernel(), good_values


def inverse_formula(target, weighted_bound):
    labels_inverse = []
    monomials = []
    for y_power in range(2):
        for x_power in range(weighted_bound // 2 + 1):
            if 2 * x_power + 3 * y_power > weighted_bound:
                continue
            labels_inverse.append((x_power, y_power))
            monomial = h_pow(X_H, x_power)
            if y_power:
                monomial = h_mul(monomial, Y_H)
            monomials.append(monomial)
    kernel, good_values = h_relation_kernel(
        tuple(monomials) + tuple(h_mul(target, value) for value in monomials)
    )
    split = len(monomials)
    for candidate in kernel.basis():
        candidate = vector(quadratic, candidate)
        numerator = h_scale(ONE_H, 0)
        denominator = h_scale(ONE_H, 0)
        for index, monomial in enumerate(monomials):
            numerator = h_add(numerator, h_scale(monomial, candidate[index]))
            denominator = h_add(
                denominator, h_scale(monomial, candidate[split + index])
            )
        if h_is_zero_mod_19_power(denominator):
            continue
        if not h_is_zero_mod_19_power(
            h_add(numerator, h_mul(target, denominator))
        ):
            continue
        pivot = next(
            (value for value in candidate if value.valuation() == 0), None
        )
        if pivot is None:
            continue
        candidate /= pivot
        return {
            "weighted_bound": weighted_bound,
            "kernel_dimension": int(kernel.dimension()),
            "good_W_residues_used": good_values,
            "monomials_X_power_Y_power": [list(label) for label in labels_inverse],
            "numerator_coefficients_mod_19_power_1_omega": [
                residue_coordinates(value) for value in candidate[:split]
            ],
            "denominator_coefficients_mod_19_power_1_omega": [
                residue_coordinates(value) for value in candidate[split:]
            ],
            "formula": "target=-numerator/denominator",
            "literal_cubic_algebra_replay_mod_19_power": True,
        }
    raise ArithmeticError(
        f"no p-adic inverse identity at weighted bound {weighted_bound}; "
        f"kernel dimension {kernel.dimension()}"
    )


inverse_W = inverse_formula(W_H, 4)
inverse_z = inverse_formula(Z_H, 10)

output = {
    "schema": "elkies-k3.q80-third-q12-riemann-roch-p19-adic-sample.v2",
    "status": "PASS_EXACT_THIRD_Q12_WEIERSTRASS_P19_ADIC_SAMPLE",
    "specialization": {
        "u": "-2",
        "prime": prime,
        "digits": digits,
        "modulus": modulus,
        "base_U_coefficients_1_omega": [base_constant, base_anti],
        "legacy_T1_reduction": is_legacy_sample,
    },
    "infinity": {
        "simple_branch_residue_coefficients_1_omega": [
            value % prime for value in residue_coordinates(simple_seed)
        ],
        "double_branch_residue_coefficients_1_omega": [
            value % prime for value in residue_coordinates(double_seed)
        ],
        "double_branch_parameter_extension_degree": branch_degree,
        "integral_basis_valuations_simple_double": integral_basis_valuations,
    },
    "riemann_roch": {
        "basis": ["1", "z", "(z^2+A*z+B)/(L*Q)"],
        "degree_bound": degree_bound,
        "labels_basis_index_W_degree": [list(label) for label in labels],
        "dimensions_L0_to_L3": dimensions,
        "x_coordinates_mod_19_power_1_omega": [
            residue_coordinates(value) for value in x_quadratic_coordinates
        ],
        "y_coordinates_mod_19_power_1_omega": [
            residue_coordinates(value) for value in y_quadratic_coordinates
        ],
        "pole_orders_x_y": [2, 3],
        "gauge": (
            "at t=1/W and the derived simple infinity-branch lift: "
            "x=t^-2+... with constant 0; "
            "y=t^-3+0*t^-2+... with constant 0"
        ),
    },
    "weierstrass": {
        "relation_coefficients_Y2_XY_Y_minusX3_minusX2_minusX_minus1_mod_19_power_1_omega": [
            residue_coordinates(value) for value in relation
        ],
        "a1_a2_a3_a4_a6_mod_19_power_1_omega": [
            residue_coordinates(value) for value in weierstrass_coefficients
        ],
        "discriminant_mod_19_power_1_omega": residue_coordinates(curve.discriminant()),
        "j_mod_19_power_1_omega": residue_coordinates(curve.j_invariant()),
        "literal_substitution_into_cubic_algebra_mod_19_power": True,
        "legacy_T1_reduction_replayed": is_legacy_sample,
    },
    "birational_maps": {
        "forward": {
            "X_weierstrass_as_function_of_W_z": function_record(weierstrass_x),
            "Y_weierstrass_as_function_of_W_z": function_record(weierstrass_y),
            "literal_substitution_into_weierstrass_mod_19_power": True,
        },
        "inverse": {
            "W": inverse_W,
            "z": inverse_z,
            "literal_cubic_algebra_replay_mod_19_power": True,
        },
    },
    "inputs": {
        "source": {"path": str(args.source.relative_to(ROOT)), "sha256": sha256(args.source)},
        "lift": {"path": str(args.lift.relative_to(ROOT)), "sha256": sha256(args.lift)},
        "basis": {"path": str(args.basis.relative_to(ROOT)), "sha256": sha256(args.basis)},
        "control": {"path": str(args.control.relative_to(ROOT)), "sha256": sha256(args.control)},
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "the p-adic infinity branches and integral-basis valuations at one exact-gauge base value",
            "dimensions 1,1,2,3 for L(0P) through L(3P) modulo 19^5",
            "normalized pole-two and pole-three generators at that p-adic sample",
            "the long Weierstrass equation and forward map modulo 19^5 at that sample",
            "inverse formulas for W and z at the finite-field pinned weighted bounds",
            "literal replay of the relation in the cubic algebra",
            "reduction to the legacy T=1 positive control when U is in its aligned residue class",
        ],
        "not_proved": [
            "generic p-adic Weierstrass relation/maps",
            "characteristic-zero reconstruction, minimization, fibres, or marking",
        ],
    },
    "reproduce": (
        "sage -python elkies-k3/scripts/compile_q80_third_q12_riemann_roch_p19_adic_sample.sage "
        f"--base-constant {base_constant} --base-anti {base_anti}"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"p-adic Riemann--Roch artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12PADICRR|p=19|digits=5|U={base_constant}+{base_anti}omega|"
    f"RR=1,1,2,3|poles=2,3|weierstrass=1|maps=both|legacy={int(is_legacy_sample)}|"
    "status=PASS_EXACT_THIRD_Q12_WEIERSTRASS_P19_ADIC_SAMPLE",
    flush=True,
)
