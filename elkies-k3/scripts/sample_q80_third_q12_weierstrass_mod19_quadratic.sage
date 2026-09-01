#!/usr/bin/env sage -python
"""Convert one resolved third-q12 fibre to Weierstrass form over GF(19^2).

The weighted infinity cubic has the rational simple branch ``xi=-6`` and a
double branch ``xi=3``.  This worker computes the finite maximal order of the
cubic function field, imposes regularity at the double branch by exact local
series, recovers ``L(2P)`` and ``L(3P)`` at the simple branch, and derives the
unique seven-term Weierstrass relation.  It is the sampling primitive for a
generic child-Jacobian interpolation.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    EllipticCurve,
    GF,
    FunctionField,
    LaurentSeriesRing,
    Matrix,
    PolynomialRing,
)


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "artifacts/generated-results/q80-third-q12-um2-p19-resolved-pencil.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-um2-p19-weierstrass-sample.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--new-base", type=int, default=1)
parser.add_argument(
    "--new-base-r",
    type=int,
    default=0,
    help="coefficient of r in the GF(19^2) new-base value",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


payload = json.loads(INPUT.read_text())
if payload["status"] != "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC":
    raise ValueError("resolved pencil is not certified")

base_finite = GF(19)
modulus_ring = PolynomialRing(base_finite, "m")
m = modulus_ring.gen()
finite = GF(19**2, "r", modulus=m**2 + 12 * m + 3)
r = finite.gen()
new_base = finite(args.new_base) + finite(args.new_base_r) * r


def field_element(coordinates):
    return finite(coordinates[0]) + finite(coordinates[1]) * r


# Cubic function field in the old x-coordinate over GF(19^2)(W).
old_function = FunctionField(finite, "W")
W = old_function.gen()
x_ring = PolynomialRing(old_function, "x")
x = x_ring.gen()
equation = x_ring.zero()
for t_degree, w_degree, x_degree, coordinates in payload["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    equation += (
        field_element(coordinates)
        * new_base**t_degree
        * W**w_degree
        * x**x_degree
    )
if equation.degree() != 3 or not equation.is_irreducible():
    raise ArithmeticError("specialized moving cubic is not irreducible")
curve_function = old_function.extension(equation.monic(), "z")
z = curve_function.gen()
finite_order = curve_function.maximal_order()
integral_basis = finite_order.basis()
if len(integral_basis) != 3:
    raise ArithmeticError("unexpected cubic integral-basis rank")

# Weighted infinity chart t=1/W, xi=t^3*x.
infinity_ring = PolynomialRing(finite, names=("t", "xi"))
t, xi = infinity_ring.gens()
infinity_equation = infinity_ring.zero()
for t_degree, w_degree, x_degree, coordinates in payload["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    infinity_equation += (
        field_element(coordinates)
        * new_base**t_degree
        * t ** (9 - w_degree - 3 * x_degree)
        * xi**x_degree
    )
infinity_cubic = infinity_equation(t=0)
expected = infinity_cubic[xi**3] * (xi + 6) * (xi + 16) ** 2
if infinity_cubic != expected:
    raise ArithmeticError("weighted infinity branches are not -6 and -16(double)")

precision = 60
simple_series_ring = LaurentSeriesRing(finite, "s", default_prec=precision)
s = simple_series_ring.gen()
simple_x_ring = PolynomialRing(simple_series_ring, "X")
Xseries = simple_x_ring.gen()
simple_equation = simple_x_ring(infinity_equation(s, Xseries))
simple_root = simple_series_ring(-6)
for unused in range(8):
    simple_root -= simple_equation(simple_root) / simple_equation.derivative()(simple_root)
if simple_equation(simple_root).valuation() < precision - 5:
    raise ArithmeticError("simple infinity branch did not converge")

# At the double branch t=s^2 and
# xi=3+c2*s^2+c3*s^3+....  Determine c2 and the quadratic equation of c3
# from the specialized curve rather than pinning the T=1 unit.
branch_parameter_ring = PolynomialRing(finite, names=("c2", "c3"))
c2_variable, c3_variable = branch_parameter_ring.gens()
trial_series_ring = PolynomialRing(branch_parameter_ring, "s_trial")
s_trial = trial_series_ring.gen()
trial = trial_series_ring(
    infinity_equation(
        s_trial**2,
        3 + c2_variable * s_trial**2 + c3_variable * s_trial**3,
    )
)
c2_ring = PolynomialRing(finite, "c2")
c2_polynomial = c2_ring(trial[4].subs({c3_variable: 0}))
c2_roots = c2_polynomial.roots(multiplicities=False)
if len(c2_roots) != 1:
    raise ArithmeticError("double branch has no unique quadratic shift")
c2_value = c2_roots[0]
c3_ring = PolynomialRing(finite, "c3")
c3_polynomial = None
for exponent in range(5, 10):
    candidate = c3_ring(trial[exponent].subs({c2_variable: c2_value}))
    if candidate:
        c3_polynomial = candidate.monic().squarefree_part()
        break
if c3_polynomial is None or c3_polynomial.degree() != 2:
    raise ArithmeticError("double branch has no quadratic cubic-shift equation")
c3_roots = c3_polynomial.roots(multiplicities=False)
if c3_roots:
    branch_field = finite
    c3_value = c3_roots[0]
    branch_degree = 1
else:
    branch_field = finite.extension(c3_polynomial, "c3")
    c3_value = branch_field.gen()
    branch_degree = 2
double_series_ring = LaurentSeriesRing(
    branch_field, "s", default_prec=2 * precision
)
sd = double_series_ring.gen()
double_polynomial_ring = PolynomialRing(branch_field, names=("t", "xi"))
td, xid = double_polynomial_ring.gens()
double_equation_polynomial = double_polynomial_ring.zero()
for exponent, coefficient in infinity_equation.dict().items():
    double_equation_polynomial += (
        branch_field(coefficient) * td**exponent[0] * xid**exponent[1]
    )
double_x_ring = PolynomialRing(double_series_ring, "X")
Xd = double_x_ring.gen()
double_equation = double_x_ring(double_equation_polynomial(sd**2, Xd))
double_root = double_series_ring(
    3 + branch_field(c2_value) * sd**2 + c3_value * sd**3
)
for unused in range(9):
    double_root -= double_equation(double_root) / double_equation.derivative()(double_root)
if double_equation(double_root).valuation() < 2 * precision - 8:
    raise ArithmeticError("double infinity branch did not converge")


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
        for index, coefficient in enumerate(curve_function(value).list())
    )


simple_basis = [
    evaluate_element(value, s**-1, simple_root * s**-3, lambda coefficient: coefficient)
    for value in integral_basis
]
double_basis = [
    evaluate_element(value, sd**-2, double_root * sd**-6, branch_field)
    for value in integral_basis
]
integral_basis_valuations = [
    [int(simple_value.valuation()), int(double_value.valuation())]
    for simple_value, double_value in zip(simple_basis, double_basis)
]
if integral_basis_valuations != [[0, 0], [-3, -6], [-2, -4]]:
    raise ArithmeticError("unexpected infinity valuations of the finite integral basis")

# The bound three is sufficient and is verified by the dimensions below.
degree_bound = 3
labels = [
    (basis_index, degree)
    for basis_index in range(3)
    for degree in range(degree_bound + 1)
]
functions = [
    W**degree * curve_function(integral_basis[basis_index])
    for basis_index, degree in labels
]
simple_columns = [
    simple_basis[basis_index] * s**-degree
    for basis_index, degree in labels
]
double_columns = [
    double_basis[basis_index] * sd ** (-2 * degree)
    for basis_index, degree in labels
]


def riemann_roch_space(pole_bound):
    rows = []
    for exponent in range(
        min(value.valuation() for value in simple_columns), -pole_bound
    ):
        rows.append([value[exponent] for value in simple_columns])
    for exponent in range(min(value.valuation() for value in double_columns), 0):
        for coordinate in range(branch_degree):
            rows.append(
                [
                    (
                        [finite(value[exponent])]
                        if branch_degree == 1
                        else list(value[exponent].list())
                    )[coordinate]
                    for value in double_columns
                ]
            )
    return Matrix(finite, rows).right_kernel()


spaces = [riemann_roch_space(bound) for bound in range(4)]
dimensions = [int(space.dimension()) for space in spaces]
if dimensions != [1, 1, 2, 3]:
    raise ArithmeticError(f"unexpected Riemann--Roch dimensions {dimensions}")


def exact_function(coordinates):
    return sum(
        (coefficient * function for coefficient, function in zip(coordinates, functions)),
        curve_function.zero(),
    )


def simple_valuation(value):
    return int(
        evaluate_element(value, s**-1, simple_root * s**-3, lambda coefficient: coefficient).valuation()
    )


x_coordinates = next(
    coordinates
    for coordinates in spaces[2].basis()
    if simple_valuation(exact_function(coordinates)) == -2
)
y_coordinates = next(
    coordinates
    for coordinates in spaces[3].basis()
    if simple_valuation(exact_function(coordinates)) == -3
)
x_function = exact_function(x_coordinates)
y_function = exact_function(y_coordinates)

# Fix the Weierstrass gauge coherently over the new-base line.  The local
# parameter is the pinned t=1/W at the rational simple branch xi=-6.
x_local = evaluate_element(
    x_function, s**-1, simple_root * s**-3, lambda coefficient: coefficient
)
x_function /= x_local[-2]
x_local = evaluate_element(
    x_function, s**-1, simple_root * s**-3, lambda coefficient: coefficient
)
x_function -= curve_function(x_local[0])
x_local = evaluate_element(
    x_function, s**-1, simple_root * s**-3, lambda coefficient: coefficient
)
if x_local[-2] != 1 or x_local[0] != 0:
    raise ArithmeticError("failed to normalize the pole-two function")

y_local = evaluate_element(
    y_function, s**-1, simple_root * s**-3, lambda coefficient: coefficient
)
y_function /= y_local[-3]
y_local = evaluate_element(
    y_function, s**-1, simple_root * s**-3, lambda coefficient: coefficient
)
y_function -= y_local[-2] * x_function
y_local = evaluate_element(
    y_function, s**-1, simple_root * s**-3, lambda coefficient: coefficient
)
y_function -= curve_function(y_local[0])
y_local = evaluate_element(
    y_function, s**-1, simple_root * s**-3, lambda coefficient: coefficient
)
if y_local[-3] != 1 or y_local[-2] != 0 or y_local[0] != 0:
    raise ArithmeticError("failed to normalize the pole-three function")

# Coefficients multiply Y^2, XY, Y, -X^3, -X^2, -X, -1.
relation_terms = (
    y_function**2,
    x_function * y_function,
    y_function,
    -x_function**3,
    -x_function**2,
    -x_function,
    -curve_function.one(),
)


def constant_relation_kernel(function_values):
    rows = []
    for component in range(3):
        values = [curve_function(value).list()[component] for value in function_values]
        polynomial_ring = values[0].denominator().parent()
        common_denominator = polynomial_ring.one()
        for value in values:
            common_denominator = common_denominator.lcm(value.denominator())
        numerators = []
        for value in values:
            cleared = value * common_denominator
            if cleared.denominator() != 1:
                raise ArithmeticError("failed to clear relation denominator")
            numerators.append(cleared.numerator())
        max_degree = max(numerator.degree() for numerator in numerators)
        for degree in range(max_degree + 1):
            rows.append([numerator[degree] for numerator in numerators])
    return Matrix(finite, rows).right_kernel()


relation_kernel = constant_relation_kernel(relation_terms)
if relation_kernel.dimension() != 1:
    raise ArithmeticError(
        f"Weierstrass relation has dimension {relation_kernel.dimension()}"
    )
relation = relation_kernel.basis()[0]
if not relation[0] or not relation[3]:
    raise ArithmeticError("degenerate Weierstrass relation")
relation /= relation[0]
cubic_scale = relation[3]
a1 = relation[1]
a2 = relation[4]
a3 = relation[2] * cubic_scale
a4 = relation[5] * cubic_scale
a6 = relation[6] * cubic_scale**2
curve = EllipticCurve(finite, [a1, a2, a3, a4, a6])
if not curve.discriminant():
    raise ArithmeticError("singular Weierstrass output")

# The output Weierstrass coordinates are X'=c*X and Y'=c*Y.
weierstrass_x = cubic_scale * x_function
weierstrass_y = cubic_scale * y_function


def inverse_formula(target, max_bound):
    """Find target=-A(X',Y')/B(X',Y') with weighted bound."""
    for bound in range(max_bound + 1):
        labels = []
        monomials = []
        for y_power in range(2):
            for x_power in range(bound // 2 + 1):
                weight = 2 * x_power + 3 * y_power
                if weight > bound:
                    continue
                labels.append((x_power, y_power))
                monomials.append(
                    weierstrass_x**x_power * weierstrass_y**y_power
                )
        kernel = constant_relation_kernel(
            tuple(monomials) + tuple(target * value for value in monomials)
        )
        for relation_candidate in kernel.basis():
            split = len(monomials)
            denominator = sum(
                (
                    relation_candidate[split + index] * value
                    for index, value in enumerate(monomials)
                ),
                curve_function.zero(),
            )
            if not denominator:
                continue
            numerator = sum(
                (
                    relation_candidate[index] * value
                    for index, value in enumerate(monomials)
                ),
                curve_function.zero(),
            )
            if target != -numerator / denominator:
                continue
            return {
                "weighted_bound": bound,
                "monomials_X_power_Y_power": [list(label) for label in labels],
                "numerator_coefficients": [
                    coordinates(value) for value in relation_candidate[:split]
                ],
                "denominator_coefficients": [
                    coordinates(value) for value in relation_candidate[split:]
                ],
                "formula": "target=-numerator/denominator",
            }
    raise ArithmeticError(f"no inverse formula through weighted bound {max_bound}")


def coordinates(value):
    result = list(finite(value).list())
    result += [base_finite.zero()] * (2 - len(result))
    return [int(result[0]), int(result[1])]


def rational_function_record(value):
    value = old_function(value)
    return {
        "numerator_coefficients_low_to_high_W": [
            coordinates(coefficient) for coefficient in value.numerator().list()
        ],
        "denominator_coefficients_low_to_high_W": [
            coordinates(coefficient) for coefficient in value.denominator().list()
        ],
    }


def function_record(value):
    return {
        "coefficients_low_to_high_old_x": [
            rational_function_record(coefficient)
            for coefficient in curve_function(value).list()
        ]
    }


inverse_W = inverse_formula(curve_function(W), 8)
inverse_x = inverse_formula(z, 12)


output = {
    "schema": "elkies-k3.q80-third-q12-weierstrass-sample-modp2.v1",
    "status": "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_MOD19_QUADRATIC",
    "specialization": {
        "u": "-2",
        "prime": 19,
        "extension_modulus": "r^2+12*r+3",
        "new_base": int(args.new_base) % 19,
        "new_base_coefficients_1_r": [
            int(args.new_base) % 19,
            int(args.new_base_r) % 19,
        ],
    },
    "infinity": {
        "simple_branch": "xi=-6",
        "double_branch": "xi=-16",
        "double_branch_parameter": (
            "t=s^2, xi=3+c2*s^2+c3*s^3+..."
        ),
        "double_branch_c2": coordinates(c2_value),
        "double_branch_c3_minimal_polynomial": str(c3_polynomial),
        "double_branch_coefficient_extension_degree": branch_degree,
        "finite_integral_basis_valuations_simple_double": integral_basis_valuations,
    },
    "riemann_roch": {
        "degree_bound": degree_bound,
        "dimensions_L0_to_L3": dimensions,
        "x_coordinates": [coordinates(value) for value in x_coordinates],
        "y_coordinates": [coordinates(value) for value in y_coordinates],
        "pole_orders_x_y": [2, 3],
        "gauge": (
            "at t=1/W,xi=-6: x=t^-2+... with constant 0; "
            "y=t^-3+0*t^-2+... with constant 0"
        ),
    },
    "weierstrass": {
        "relation": [coordinates(value) for value in relation],
        "a1_a2_a3_a4_a6": [coordinates(value) for value in (a1, a2, a3, a4, a6)],
        "discriminant": coordinates(curve.discriminant()),
        "j": coordinates(curve.j_invariant()),
        "point_count": int(curve.cardinality()),
        "frobenius_trace": int(19**2 + 1 - curve.cardinality()),
    },
    "birational_maps": {
        "forward": {
            "X_weierstrass_as_function_of_W_old_x": function_record(weierstrass_x),
            "Y_weierstrass_as_function_of_W_old_x": function_record(weierstrass_y),
            "literal_substitution_into_weierstrass": True,
        },
        "inverse": {
            "W": inverse_W,
            "old_x": inverse_x,
            "literal_function_field_replay": True,
        },
    },
    "input": {"path": str(INPUT.relative_to(ROOT)), "sha256": sha256(INPUT)},
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": (
        "This exactly converts one smooth resolved pencil member to its "
        "elliptic Jacobian with forward and inverse function-field maps. It "
        "does not yet interpolate the generic child, "
        "classify its fibres, or lift across primes/characteristic zero."
    ),
    "reproduce": (
        "sage -python "
        "elkies-k3/scripts/sample_q80_third_q12_weierstrass_mod19_quadratic.sage "
        f"--new-base {args.new_base} --new-base-r {args.new_base_r}"
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12WEIERSTRASSSAMPLE|T={new_base}|"
    "RR=1,1,2,3|poles=2,3|"
    f"Delta={curve.discriminant()}|j={curve.j_invariant()}|"
    f"points={curve.cardinality()}|"
    "status=PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_MOD19_QUADRATIC",
    flush=True,
)
