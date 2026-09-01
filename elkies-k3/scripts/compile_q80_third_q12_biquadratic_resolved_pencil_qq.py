#!/usr/bin/env python3
"""Adapt the connected third-q12 compiler to the exact biquadratic horizontal."""

import argparse
import contextlib
import hashlib
import io
import json
import re
import sys
import tempfile
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
CORE = ROOT / "elkies-k3/scripts/compile_q80_third_q12_um2_p19_resolved_pencil.sage"
LATTICE = RESULTS / "q80-d7d5-mw5-height-lattice.json"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--surface", type=Path,
    default=RESULTS / "q80-fixed-u-minus2-p19-height-shell-with-po1.json",
)
parser.add_argument(
    "--horizontal", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-horizontal-qq.json",
)
parser.add_argument(
    "--operands", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json",
)
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json",
)
args = parser.parse_args()
for name in ("surface", "horizontal", "operands", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rational(record):
    return QQ(ZZ(record["numerator"])) / ZZ(record["denominator"])


horizontal = json.loads(args.horizontal.read_text())
operands = json.loads(args.operands.read_text())
if horizontal.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_HORIZONTAL_HELDOUT_P71":
    raise ValueError("exact biquadratic horizontal is not certified")
if operands.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL":
    raise ValueError("exact biquadratic operands are not certified")

q1 = read_rational(operands["biquadratic_field"]["q1"])
q2 = read_rational(operands["biquadratic_field"]["q2"])
z_ring = PolynomialRing(QQ, "z")
z = z_ring.gen()
# A relative presentation makes PARI compute an enormous integral basis even
# with ``check=False``.  The primitive element theta=a+b has this quartic
# polynomial and gives the same field immediately in the power basis.
primitive_polynomial = z**4 - 2*(q1 + q2)*z**2 + (q1 - q2)**2
field = z_ring.quotient(primitive_polynomial, "theta")
theta = field.gen()
a = (theta**3 - (3*q1 + q2)*theta) / (2*(q2 - q1))
b = theta - a
if primitive_polynomial.degree() != 4 or a**2 != field(q1) or b**2 != field(q2):
    raise ArithmeticError("operand square classes do not generate a degree-four field")
# Sparse univariate polynomials use Sage's generic Euclidean arithmetic.  The
# dense number-field specialization delegates gcds to PARI and consequently
# requests an integral basis of this enormous presentation.
base = PolynomialRing(field, "W", sparse=True)


def component_polynomials(record):
    numerator = base([field(QQ(value)) for value in record["numerator_coefficients_low_to_high"]])
    denominator = base([field(QQ(value)) for value in record["denominator_coefficients_low_to_high"]])
    return numerator, denominator


x_records = horizontal["horizontal"]["x_basis_coefficients"]
y_records = horizontal["horizontal"]["y_basis_coefficients"]
x0_numerator, x_denominator = component_polynomials(x_records[0])
x3_numerator, x3_denominator = component_polynomials(x_records[3])
y1_numerator, y_denominator = component_polynomials(y_records[1])
y2_numerator, y2_denominator = component_polynomials(y_records[2])
if x_denominator != x3_denominator or y_denominator != y2_denominator:
    raise ArithmeticError("exact horizontal components do not have common denominators")
# Combine at the polynomial level.  Constructing separate fraction-field
# elements first asks PARI for expensive gcds of enormous number-field
# polynomials even though the denominators are already identical.
x_numerator = x0_numerator + field(a*b)*x3_numerator
y_numerator = field(a)*y1_numerator + field(b)*y2_numerator


def euclidean_gcd(left, right):
    """Polynomial gcd without the number-field factorization backend."""
    left, right = base(left), base(right)
    while right:
        unused_quotient, remainder = left.quo_rem(right)
        left, right = right, remainder
    return left / left.leading_coefficient()


def cancel_and_monic(numerator, denominator):
    common = euclidean_gcd(numerator, denominator)
    numerator, numerator_remainder = numerator.quo_rem(common)
    denominator, denominator_remainder = denominator.quo_rem(common)
    if numerator_remainder or denominator_remainder:
        raise ArithmeticError("horizontal gcd cancellation failed")
    leading = denominator.leading_coefficient()
    return numerator / leading, denominator / leading


x_numerator, x_denominator = cancel_and_monic(x_numerator, x_denominator)


def monic_polynomial_square_root(value):
    if value.degree() % 2 or value.leading_coefficient() != field.one():
        raise ArithmeticError("expected a monic even-degree square")
    root_degree = value.degree() // 2
    coefficients = [field.zero()] * (root_degree + 1)
    coefficients[root_degree] = field.one()
    for target_degree in range(2*root_degree - 1, root_degree - 1, -1):
        index = target_degree - root_degree
        known = sum(
            coefficients[left] * coefficients[target_degree - left]
            for left in range(root_degree + 1)
            if 0 <= target_degree - left <= root_degree
            and left != index and target_degree - left != index
        )
        coefficients[index] = (value[target_degree] - known) / 2
    root = base(coefficients)
    if root**2 != value:
        raise ArithmeticError("expected polynomial is not a square")
    return root


h_exact = monic_polynomial_square_root(x_denominator)
# The normalized section must have denominator h^3 in y.  Dividing the raw
# common denominator by that known target avoids a second algebraic gcd.
y_leading = y_denominator.leading_coefficient()
y_numerator /= y_leading
y_denominator /= y_leading
y_common, y_denominator_remainder = y_denominator.quo_rem(h_exact**3)
if y_denominator_remainder:
    raise ArithmeticError("raw y denominator is not divisible by h^3")
y_numerator, y_numerator_remainder = y_numerator.quo_rem(y_common)
if y_numerator_remainder:
    raise ArithmeticError("raw y numerator does not contain its expected fixed factor")
y_denominator = h_exact**3
if (x_numerator.degree(), x_denominator.degree()) != (8, 4):
    raise ArithmeticError("normalized exact horizontal x has the wrong degrees")
if (y_numerator.degree(), y_denominator.degree()) != (12, 6):
    raise ArithmeticError("normalized exact horizontal y has the wrong degrees")


def coefficient_strings(polynomial):
    return [str(coefficient) for coefficient in polynomial.list()]


normalized_horizontal = {
    "schema": "elkies-k3.q80-po0-rur-third-q12-modp.v1",
    "status": "PASS_EXACT_MODP2_THIRD_Q12_HORIZONTAL_FROBENIUS_ORBIT",
    "third_q12": {
        "candidates_up_to_sign": [
            {
                "x": {
                    "numerator_coefficients_low_to_high": coefficient_strings(x_numerator),
                    "denominator_coefficients_low_to_high": coefficient_strings(x_denominator),
                },
                "y": {
                    "numerator_coefficients_low_to_high": coefficient_strings(y_numerator),
                    "denominator_coefficients_low_to_high": coefficient_strings(y_denominator),
                },
            }
        ]
    },
}

local_root = ROOT / "artifacts/local/elkies-k3"
local_root.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(dir=local_root) as temporary_directory:
    temporary = Path(temporary_directory)
    normalized_path = temporary / "horizontal.json"
    core_output = temporary / "resolved.json"
    normalized_path.write_text(json.dumps(normalized_horizontal, indent=2, sort_keys=True) + "\n")

    source = CORE.read_text()
    path_pattern = re.compile(
        r"SURFACE = ROOT / \(.*?DEFAULT_OUTPUT = ROOT / \(.*?\n\)",
        re.DOTALL,
    )
    replacement = (
        f"SURFACE = Path({str(args.surface)!r})\n"
        f"HORIZONTAL = Path({str(normalized_path)!r})\n"
        f"LATTICE = Path({str(LATTICE)!r})\n"
        f"DEFAULT_OUTPUT = Path({str(core_output)!r})"
    )
    source, path_replacements = path_pattern.subn(replacement, source, count=1)
    field_pattern = re.compile(
        r"prime = 19\nbase_finite = GF\(prime\).*?base_ring = PolynomialRing\(finite, \"W\"\)\nW = base_ring.gen\(\)",
        re.DOTALL,
    )
    q1_text = repr(q1)
    q2_text = repr(q2)
    field_replacement = f'''prime = 0
base_finite = QQ
z_ring_exact = PolynomialRing(QQ, "z_exact")
z_exact = z_ring_exact.gen()
q1_exact = QQ({q1_text!r})
q2_exact = QQ({q2_text!r})
primitive_exact = z_exact**4 - 2*(q1_exact + q2_exact)*z_exact**2 + (q1_exact - q2_exact)**2
finite = z_ring_exact.quotient(primitive_exact, "theta")
theta = finite.gen()
a = (theta**3 - (3*q1_exact + q2_exact)*theta) / (2*(q2_exact - q1_exact))
b = theta - a
base_ring = PolynomialRing(finite, "W", sparse=True)
W = base_ring.gen()'''
    source, field_replacements = field_pattern.subn(field_replacement, source, count=1)
    coordinate_pattern = re.compile(
        r"def field_coordinates\(value\):\n.*?return \[int\(values\[0\]\), int\(values\[1\]\)\]",
        re.DOTALL,
    )
    coordinate_replacement = '''def field_coordinates(value):
    try:
        converted_value = finite(value)
    except Exception as error:
        raise TypeError(
            f"cannot encode field value type={type(value)!r} repr={value!r}"
        ) from error
    return [str(converted_value)]'''
    source, coordinate_replacements = coordinate_pattern.subn(coordinate_replacement, source, count=1)
    pole_old = '''h_factors = tuple(x_denominator.factor())
if len(h_factors) != 1 or int(h_factors[0][1]) != 2:
    raise ArithmeticError("horizontal x denominator is not one square")
h = h_factors[0][0].monic()'''
    pole_new = '''def exact_polynomial_square_root(value):
    if value.degree() % 2:
        raise ArithmeticError("horizontal x denominator has odd degree")
    root_degree = value.degree() // 2
    if value.leading_coefficient() != finite.one():
        raise ArithmeticError("horizontal x denominator is not monic")
    coefficients = [finite.zero()] * (root_degree + 1)
    coefficients[root_degree] = finite.one()
    for target_degree in range(2*root_degree - 1, root_degree - 1, -1):
        index = target_degree - root_degree
        known = sum(
            coefficients[left] * coefficients[target_degree - left]
            for left in range(root_degree + 1)
            if 0 <= target_degree - left <= root_degree
            and left != index and target_degree - left != index
        )
        coefficients[index] = (
            value[target_degree] - known
        ) / (2*coefficients[root_degree])
    root = base_ring(coefficients)
    if root**2 != value:
        raise ArithmeticError("horizontal x denominator is not a square")
    return root

h = exact_polynomial_square_root(x_denominator)
if h.degree() != 2:
    raise ArithmeticError("horizontal pole divisor does not have degree two")'''
    if (
        path_replacements != 1
        or field_replacements != 1
        or coordinate_replacements != 1
        or source.count(pole_old) != 1
    ):
        raise ArithmeticError("immutable core no longer matches the exact adapter contract")
    source = source.replace(pole_old, pole_new, 1)
    smith_pattern = re.compile(
        r"smith, smith_left, smith_right = module_matrix\.smith_form\(\)\n"
        r".*?saturated_module = smith_left\.inverse\(\)\[:, :3\]",
        re.DOTALL,
    )
    smith_replacement = '''def residue_constraint_matrix(basis):
    residue_degree = h.degree()
    rows = []
    for module_row in range(basis.nrows()):
        for output_degree in range(residue_degree):
            row = []
            for column in range(basis.ncols()):
                for coefficient_degree in range(residue_degree):
                    residue = (basis[module_row, column] * W**coefficient_degree).mod(h)
                    row.append(residue[output_degree])
            rows.append(row)
    return Matrix(finite, rows)


def saturate_once_at_h(basis):
    residue_degree = h.degree()
    constraint = residue_constraint_matrix(basis)
    kernel = constraint.right_kernel()
    if constraint.rank() != 4 or kernel.dimension() != 2:
        raise ArithmeticError("unexpected h-saturation residue kernel")
    chosen = None
    for kernel_vector in kernel.basis():
        coefficients = [
            base_ring(list(kernel_vector[
                column*residue_degree:(column + 1)*residue_degree
            ]))
            for column in range(basis.ncols())
        ]
        pivot = next(
            (
                column for column, coefficient in enumerate(coefficients)
                if coefficient and coefficient.degree() == 0
            ),
            None,
        )
        if pivot is not None:
            chosen = coefficients, pivot
            break
    if chosen is None:
        raise ArithmeticError("h-saturation kernel has no constant pivot")
    coefficients, pivot = chosen
    scale = coefficients[pivot][0]
    coefficients = [coefficient / scale for coefficient in coefficients]
    combined = sum(
        (
            coefficient * basis.column(column)
            for column, coefficient in enumerate(coefficients)
        ),
        vector(base_ring, [0] * basis.nrows()),
    )
    quotient_entries = []
    for entry in combined:
        quotient, remainder = entry.quo_rem(h)
        if remainder:
            raise ArithmeticError("h-saturation vector is not divisible by h")
        quotient_entries.append(quotient)
    enlarged = Matrix(base_ring, basis)
    enlarged.set_column(pivot, vector(base_ring, quotient_entries))
    return enlarged


# Away from h the rows (1,2,4) already have determinant h^7, so h is the
# only possible determinantal support.  The original divisor is h^3; three
# checked index-h enlargements remove it completely.
saturated_module = Matrix(base_ring, module_matrix)
for unused_saturation_step in range(3):
    saturated_module = saturate_once_at_h(saturated_module)
final_residue_constraint = residue_constraint_matrix(saturated_module)
if final_residue_constraint.rank() != 6 or final_residue_constraint.right_kernel().dimension() != 0:
    raise ArithmeticError("three h-saturations did not produce a primitive module")
smith_degrees = (0, 0, 6)'''
    source, smith_replacements = smith_pattern.subn(smith_replacement, source, count=1)
    if smith_replacements != 1:
        raise ArithmeticError("immutable core Smith block no longer matches adapter contract")
    source = source.replace(
        'raise ArithmeticError("unexpected shifted Popov weights")',
        'raise ArithmeticError(f"unexpected shifted Popov weights: {generator_weights}; pivots={popov_pivots}")',
        1,
    )
    popov_line = "popov_module = shifted_weak_popov(saturated_module)"
    popov_canonical = '''popov_module = shifted_weak_popov(saturated_module)
provisional_pivots = tuple(
    column_pivot(popov_module.column(index))
    for index in range(popov_module.ncols())
)
provisional_weights = tuple(degree - 8 for degree, unused in provisional_pivots)
if sorted(provisional_weights) != [0, 0, 2]:
    raise ArithmeticError(f"unexpected shifted Popov weight multiset: {provisional_weights}")
positive_column = provisional_weights.index(2)
zero_columns = [index for index, weight in enumerate(provisional_weights) if weight == 0]
canonical_order = [zero_columns[0], positive_column, zero_columns[1]]
popov_module = popov_module.matrix_from_columns(canonical_order)'''
    if source.count(popov_line) != 1:
        raise ArithmeticError("immutable core Popov call no longer matches adapter contract")
    source = source.replace(popov_line, popov_canonical, 1)
    fibre_pattern = re.compile(
        r"delta = 4 \* A\*\*3 \+ 27 \* B\*\*2\n"
        r"star_factor = next\(factor for factor, exponent in delta\.factor\(\) if exponent == 7\)\n"
        r"star_root = -star_factor\[0\] / star_factor\[1\]\n"
        r"node_ring = PolynomialRing\(finite, \"x_node\"\)\n"
        r"x_node = node_ring\.gen\(\)\n"
        r"node_cubic = x_node\*\*3 \+ A\(star_root\) \* x_node \+ B\(star_root\)\n"
        r"singular_roots = node_cubic\.gcd\(node_cubic\.derivative\(\)\)\.roots\(\n"
        r"    multiplicities=False\n"
        r"\)\n"
        r"if len\(singular_roots\) != 1:\n"
        r"    raise ArithmeticError\(\"finite I1\* cubic has no unique singular x\"\)\n"
        r"singular_x = singular_roots\[0\]",
    )
    fibre_replacement = '''rational_base_ring = PolynomialRing(QQ, "W_rational")
def quotient_constant_to_rational(coefficient):
    lifted = finite(coefficient).lift()
    if lifted.degree() > 0:
        raise ArithmeticError("parent coefficient is not rational")
    return QQ(lifted[0])

rational_A = rational_base_ring([quotient_constant_to_rational(coefficient) for coefficient in A.list()])
rational_B = rational_base_ring([quotient_constant_to_rational(coefficient) for coefficient in B.list()])
rational_delta = 4*rational_A**3 + 27*rational_B**2
rational_star_factor = next(
    factor for factor, exponent in rational_delta.factor() if exponent == 7
)
rational_star_root = -rational_star_factor[0] / rational_star_factor[1]
star_root = finite(rational_star_root)
rational_node_ring = PolynomialRing(QQ, "x_node_rational")
x_node_rational = rational_node_ring.gen()
rational_node_cubic = (
    x_node_rational**3
    + rational_A(rational_star_root)*x_node_rational
    + rational_B(rational_star_root)
)
singular_factor = rational_node_cubic.gcd(rational_node_cubic.derivative()).monic()
if singular_factor.degree() == 1:
    rational_singular_x = -singular_factor[0] / singular_factor[1]
elif singular_factor.degree() == 2 and singular_factor.discriminant() == 0:
    rational_singular_x = -singular_factor[1] / (2*singular_factor[2])
else:
    raise ArithmeticError("finite I1* cubic has no unique singular x")
if singular_factor(rational_singular_x):
    raise ArithmeticError("finite I1* singular root replay failed")
singular_x = finite(rational_singular_x)'''
    source, fibre_replacements = fibre_pattern.subn(fibre_replacement, source, count=1)
    if fibre_replacements != 1:
        raise ArithmeticError("immutable core D5 factor block no longer matches adapter contract")
    elimination_pattern = re.compile(
        r"nonzero_coefficients = tuple\(value for value in eliminated\.list\(\) if value\)\n"
        r"fixed_factor = nonzero_coefficients\[0\]\n"
        r"for coefficient in nonzero_coefficients\[1:\]:\n"
        r"    fixed_factor = fixed_factor\.gcd\(coefficient\)\n"
        r"expected_fixed = plane_ring\(h\(W0\) \*\* 2 \* x0 - Nx\(W0\)\)\n"
        r"if fixed_factor \* expected_fixed\.leading_coefficient\(\) != \(\n"
        r"    expected_fixed \* fixed_factor\.leading_coefficient\(\)\n"
        r"\):\n"
        r"    raise ArithmeticError\(\"elimination fixed factor is not h\^2\*x-Nx\"\)\n"
        r"moving = parameter_polynomials\(\n"
        r"    \[coefficient // fixed_factor for coefficient in eliminated\.list\(\)\]\n"
        r"\)",
    )
    elimination_replacement = '''expected_fixed = plane_ring(h(W0) ** 2 * x0 - Nx(W0))
def divide_by_horizontal_factor(coefficient):
    x_degree = coefficient.degree(x0)
    coefficients = [base_ring.zero()] * (x_degree + 1)
    for exponent, scalar in coefficient.dict().items():
        w_degree, current_x_degree = exponent
        coefficients[current_x_degree] += finite(scalar) * W**w_degree
    quotient_coefficients = [base_ring.zero()] * x_degree
    divisor_leading = h**2
    divisor_constant = -Nx
    for current_x_degree in range(x_degree, 0, -1):
        quotient, remainder = coefficients[current_x_degree].quo_rem(divisor_leading)
        if remainder:
            raise ArithmeticError("elimination leading coefficient is not divisible by h^2")
        quotient_coefficients[current_x_degree - 1] = quotient
        coefficients[current_x_degree] = base_ring.zero()
        coefficients[current_x_degree - 1] -= quotient * divisor_constant
    if coefficients[0]:
        raise ArithmeticError("elimination has a nonzero horizontal-factor remainder")
    answer = sum(
        plane_ring(quotient(W0)) * x0**degree
        for degree, quotient in enumerate(quotient_coefficients)
    )
    if answer * expected_fixed != coefficient:
        raise ArithmeticError("horizontal-factor synthetic division replay failed")
    return answer

moving_coefficients = []
for coefficient in eliminated.list():
    moving_coefficients.append(divide_by_horizontal_factor(coefficient))
fixed_factor = expected_fixed
moving = parameter_polynomials(moving_coefficients)'''
    source, elimination_replacements = elimination_pattern.subn(
        elimination_replacement, source, count=1
    )
    if elimination_replacements != 1:
        raise ArithmeticError("immutable core elimination block no longer matches adapter contract")
    residue_lines = '''    numerator = local_ring(value.numerator())
    denominator = local_ring(value.denominator())'''
    residue_normalization = '''    numerator = local_ring(value.numerator())
    denominator = local_ring(value.denominator())
    # Fraction reduction over an explicit algebraic quotient field does not
    # always cancel a common exceptional monomial.  Remove the common t-power
    # literally before testing whether the remaining denominator is a unit.
    numerator_t_order = min(exponent[0] for exponent in numerator.dict())
    denominator_t_order = min(exponent[0] for exponent in denominator.dict())
    common_t_order = min(numerator_t_order, denominator_t_order)
    if common_t_order:
        numerator = local_ring({
            (exponent[0] - common_t_order, exponent[1], exponent[2]): coefficient
            for exponent, coefficient in numerator.dict().items()
        })
        denominator = local_ring({
            (exponent[0] - common_t_order, exponent[1], exponent[2]): coefficient
            for exponent, coefficient in denominator.dict().items()
        })'''
    if source.count(residue_lines) != 1:
        raise ArithmeticError("immutable core residue normalization no longer matches adapter contract")
    source = source.replace(residue_lines, residue_normalization, 1)
    saved_argv = sys.argv
    namespace = {"__file__": str(CORE), "__name__": "__main__"}
    sys.argv = [str(CORE), "--output", str(core_output)]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(source, str(CORE), "exec"), namespace)
    finally:
        sys.argv = saved_argv
    compiled = json.loads(core_output.read_text())

if compiled.get("status") != "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC":
    raise ArithmeticError("adapted core did not pass its exact connected-quotient gates")
compiled["schema"] = "elkies-k3.q80-third-q12-biquadratic-resolved-pencil-qq.v1"
compiled["status"] = "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_RESOLVED_PENCIL"
compiled["specialization"] = {
    "u": "-2",
    "field_degree": 4,
    "primitive_element": "theta=a+b",
    "primitive_polynomial": "theta^4-2*(q1+q2)*theta^2+(q1-q2)^2",
    "coefficient_encoding": "one exact Sage expression in theta per coefficient",
}
compiled["inputs"] = {
    "surface": {"path": str(args.surface.relative_to(ROOT)), "sha256": sha256(args.surface)},
    "horizontal": {"path": str(args.horizontal.relative_to(ROOT)), "sha256": sha256(args.horizontal)},
    "operands": {"path": str(args.operands.relative_to(ROOT)), "sha256": sha256(args.operands)},
    "lattice": {"path": str(LATTICE.relative_to(ROOT)), "sha256": sha256(LATTICE)},
    "core": {"path": str(CORE.relative_to(ROOT)), "sha256": sha256(CORE)},
    "adapter": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
}
compiled["claim_boundary"] = {
    "proved": [
        "exact biquadratic horizontal literal replay on the characteristic-zero parent",
        "complete characteristic-zero Smith saturation and shifted-Popov ambient",
        "resolved D7 complete ideal and finite D5 connected quotient",
        "rank-five gate with a two-dimensional pencil",
        "exact moving equation of degrees (2,9,3) over the biquadratic field",
    ],
    "not_proved": [
        "generic genus one over characteristic zero",
        "a minimal characteristic-zero child Jacobian or birational maps",
        "the characteristic-zero A5+A3+3A1 marking or Mordell--Weil rank",
    ],
}
compiled["reproduce"] = (
    "sage -python elkies-k3/scripts/compile_q80_third_q12_biquadratic_resolved_pencil_qq.py "
    f"--output {args.output}"
)
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(compiled, indent=2, sort_keys=True) + "\n")
print(
    "Q80THIRDQ12BIQUADRATICPENCIL|field_degree=4|smith=0,0,6|ambient=7|"
    "gate_rank=5|kernel=2|moving_degrees=2,9,3|"
    "status=PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_RESOLVED_PENCIL",
    flush=True,
)
