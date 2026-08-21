#!/usr/bin/env sage
"""Test whether the modular q80 surface parameter has a rigid QQ Hensel lift.

This is a structural probe for the missing characteristic-zero algebraization
of the slope-8/87 branch.  The bounded GF(7) normalization gives rational
surface functions of degrees ``5/4,10/8,8/6,15/12``.  After fixing the two
PGL2 parameter freedoms at the CM point, this script asks whether matching
the four local surface series gives a full-rank coefficient system modulo 7.

With ``--qq-series``, the same rigid system is lifted 7-adically and rationally
reconstructed, then checked against every supplied characteristic-zero jet
coefficient.  The resulting rational functions remain candidates for a global
parameterization until a separate direct substitution verifies them.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing, QQ, ZZ, Zmod, vector


PARAMETER_PATH = Path(
    "artifacts/generated-results/q80-cm24-slope-8-87-gf7-parameter.json"
)
parser = argparse.ArgumentParser()
parser.add_argument("--order", type=int, default=30)
parser.add_argument(
    "--qq-series",
    help="exact CRT-reconstructed QQ surface-series artifact to Hensel lift",
)
parser.add_argument(
    "--qq-output",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-parameter.json",
)
parser.add_argument("--max-hensel-digits", type=int, default=260)
arguments = parser.parse_args()
ORDER = arguments.order
if ORDER < 4:
    parser.error("--order must be at least 4")
DEGREES = {
    "D": (5, 4),
    "P": (10, 8),
    "Q": (8, 6),
    "E": (15, 12),
}

finite = GF(7)
polynomials = PolynomialRing(finite, "t")
t = polynomials.gen()
functions = polynomials.fraction_field()
artifact = json.loads(PARAMETER_PATH.read_text())
assert artifact["field"] == "GF(7)"


def truncate_product(left, right, order=ORDER):
    result = [left[0].parent()(0) for _ in range(order)]
    for i, left_value in enumerate(left[:order]):
        for j, right_value in enumerate(right[:order-i]):
            result[i+j] += left_value*right_value
    return result


def invert_series(series, order=ORDER):
    assert series[0]
    result = [series[0]**-1]+[series[0].parent()(0)]*(order-1)
    for degree in range(1, order):
        result[degree] = -sum(
            series[index]*result[degree-index]
            for index in range(1, min(degree, len(series)-1)+1)
        )/series[0]
    return result


def rational_series(numerator, denominator, order=ORDER):
    numerator_series = list(numerator)+[numerator[0].parent()(0)]*order
    denominator_series = list(denominator)+[denominator[0].parent()(0)]*order
    return truncate_product(
        numerator_series[:order], invert_series(denominator_series, order), order
    )


def compose_series(outer, inner, order=ORDER):
    result = [outer[0].parent()(0)]*order
    power = [outer[0].parent()(1)]+[outer[0].parent()(0)]*(order-1)
    for coefficient in outer[:order]:
        for index, value in enumerate(power):
            result[index] += coefficient*value
        power = truncate_product(power, inner, order)
    return result


def normalized_coefficients(value, degrees):
    numerator_degree, denominator_degree = degrees
    value = functions(value)
    numerator = polynomials(value.numerator())
    denominator = polynomials(value.denominator())
    assert numerator.degree() <= numerator_degree
    assert denominator.degree() <= denominator_degree
    scale = denominator[0]**-1
    numerator *= scale
    denominator *= scale
    return (
        [numerator[index] for index in range(1, numerator_degree+1)],
        [denominator[index] for index in range(1, denominator_degree+1)],
    )


# Keep the osculating parameter supplied by the normalization artifact.  Its
# local P function has normalized numerator coefficient -1 and denominator
# coefficient -2.  Fixing those two coefficients removes the PGL2 freedom
# without increasing any of the certified numerator/denominator degrees.
old_parameter = functions(t)
parameter_functions = {}
seed_blocks = {}
for name, degrees in DEGREES.items():
    value = functions(artifact["functions"][name]["value"])(old_parameter)
    parameter_functions[name] = value
    seed_blocks[name] = normalized_coefficients(value, degrees)


def block_series(block, degrees, order=ORDER):
    numerator_degree, denominator_degree = degrees
    numerator_coefficients, denominator_coefficients = block
    numerator = [finite(0)]+list(numerator_coefficients)
    denominator = [finite(1)]+list(denominator_coefficients)
    assert len(numerator) == numerator_degree+1
    assert len(denominator) == denominator_degree+1
    return rational_series(numerator, denominator, order)


seed_parameter_series = block_series(seed_blocks["P"], DEGREES["P"])
assert seed_parameter_series[1] == -1
assert seed_blocks["P"][1][0] == -2


def revert_series(series, order=ORDER):
    assert series[0] == 0 and series[1]
    inverse = [finite(0)]*order
    inverse[1] = series[1]**-1
    for degree in range(2, order):
        trial = compose_series(series, inverse, degree+1)[degree]
        inverse[degree] = -trial/series[1]
    identity = compose_series(series, inverse, order)
    assert identity[0] == 0 and identity[1] == 1
    assert not any(identity[2:])
    return inverse


parameter_inverse = revert_series(seed_parameter_series)
branch_series = {"P": [finite(0), finite(1)]+[finite(0)]*(ORDER-2)}
for name in ("D", "Q", "E"):
    parameter_value = block_series(seed_blocks[name], DEGREES[name])
    branch_series[name] = compose_series(
        parameter_value, parameter_inverse, ORDER
    )


block_order = ("D", "P", "Q", "E")
offsets = {}
seed = []
for name in block_order:
    numerator, denominator = seed_blocks[name]
    start = len(seed)
    seed.extend(numerator)
    seed.extend(denominator)
    offsets[name] = (start, len(numerator), len(denominator))


def unpack(values, name):
    start, numerator_length, denominator_length = offsets[name]
    numerator = values[start:start+numerator_length]
    denominator = values[
        start+numerator_length:start+numerator_length+denominator_length
    ]
    return numerator, denominator


def equations(values):
    blocks = {name: unpack(values, name) for name in block_order}
    parameter_series = block_series(blocks["P"], DEGREES["P"])
    result = [blocks["P"][0][0]+1, blocks["P"][1][0]+2]
    for name in ("D", "Q", "E"):
        numerator, denominator = blocks[name]
        numerator_series = [finite(0)]+list(numerator)+[finite(0)]*ORDER
        denominator_series = [finite(1)]+list(denominator)+[finite(0)]*ORDER
        pulled_back = compose_series(
            branch_series[name], parameter_series, ORDER
        )
        residual = truncate_product(
            denominator_series[:ORDER], pulled_back, ORDER
        )
        residual = [
            value-numerator_series[index]
            for index, value in enumerate(residual)
        ]
        result.extend(residual[1:])
    return vector(finite, result)


seed = vector(finite, seed)
seed_residual = equations(seed)
assert not any(seed_residual)


class Dual:
    def __init__(self, value, derivative=0):
        self.value = finite(value)
        self.derivative = finite(derivative)

    def parent(self):
        return finite

    def __add__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return Dual(self.value+other.value, self.derivative+other.derivative)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, -self.derivative)

    def __sub__(self, other):
        return self+(-other if isinstance(other, Dual) else Dual(-other))

    def __rsub__(self, other):
        return Dual(other)-self

    def __mul__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return Dual(
            self.value*other.value,
            self.derivative*other.value+self.value*other.derivative,
        )

    __rmul__ = __mul__

    def __pow__(self, exponent):
        if exponent == -1:
            return Dual(
                self.value**-1,
                -self.derivative/(self.value**2),
            )
        if exponent == 0:
            return Dual(1)
        result = Dual(1)
        for _ in range(exponent):
            result *= self
        return result

    def __truediv__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return self*(other**-1)

    def __rtruediv__(self, other):
        return Dual(other)/self


def dual_equations(values):
    # Reuse the same arithmetic with dual-valued inputs.  The fixed branch
    # coefficients are promoted by the overloaded operators as needed.
    blocks = {name: unpack(values, name) for name in block_order}
    parameter_series = block_series(blocks["P"], DEGREES["P"])
    result = [blocks["P"][0][0]+1, blocks["P"][1][0]+2]
    for name in ("D", "Q", "E"):
        numerator, denominator = blocks[name]
        numerator_series = [Dual(0)]+list(numerator)+[Dual(0)]*ORDER
        denominator_series = [Dual(1)]+list(denominator)+[Dual(0)]*ORDER
        promoted_branch = [Dual(value) for value in branch_series[name]]
        pulled_back = compose_series(
            promoted_branch, parameter_series, ORDER
        )
        residual = truncate_product(
            denominator_series[:ORDER], pulled_back, ORDER
        )
        residual = [
            value-numerator_series[index]
            for index, value in enumerate(residual)
        ]
        result.extend(residual[1:])
    return result


jacobian_columns = []
for column in range(len(seed)):
    dual_seed = [Dual(value, index == column) for index, value in enumerate(seed)]
    jacobian_columns.append(
        vector(finite, [value.derivative for value in dual_equations(dual_seed)])
    )
jacobian = Matrix(finite, jacobian_columns).transpose()

status = "PASS_RIGID_MOD7" if jacobian.rank() == len(seed) else "FAIL_NONRIGID_MOD7"
print(
    "Q80SURFACEPARAM|field=GF(7)|slope=8/87|"
    f"order={ORDER}|unknowns={len(seed)}|equations={len(seed_residual)}|"
    f"jacobian_rank={jacobian.rank()}|nullity={len(seed)-jacobian.rank()}|"
    f"gauge=P_num_t:-1,P_den_t:-2|status={status}",
    flush=True,
)

assert jacobian.rank() == len(seed)


def equations_over(values, exact_branch, base):
    """Evaluate the coefficient-matching system over QQ or Z/7^n."""
    blocks = {name: unpack(values, name) for name in block_order}

    def series_for_block(block, degrees):
        numerator_degree, denominator_degree = degrees
        numerator_values, denominator_values = block
        numerator = [base(0)]+list(numerator_values)
        denominator = [base(1)]+list(denominator_values)
        assert len(numerator) == numerator_degree+1
        assert len(denominator) == denominator_degree+1
        return rational_series(numerator, denominator, ORDER)

    parameter_series = series_for_block(blocks["P"], DEGREES["P"])
    result = [blocks["P"][0][0]+base(1), blocks["P"][1][0]+base(2)]
    for name in ("D", "Q", "E"):
        numerator, denominator = blocks[name]
        numerator_series = [base(0)]+list(numerator)+[base(0)]*ORDER
        denominator_series = [base(1)]+list(denominator)+[base(0)]*ORDER
        pulled_back = compose_series(
            exact_branch[name], parameter_series, ORDER
        )
        residual = truncate_product(
            denominator_series[:ORDER], pulled_back, ORDER
        )
        result.extend(
            residual[index]-numerator_series[index]
            for index in range(1, ORDER)
        )
    return result


if arguments.qq_series:
    series_path = Path(arguments.qq_series)
    exact_artifact = json.loads(series_path.read_text())
    if exact_artifact.get("schema") != "q80-cm24-qq-surface-series-v1":
        raise ValueError("unexpected exact-series schema")
    if exact_artifact.get("slope") != "8/87":
        raise ValueError("exact series belongs to another branch")
    if exact_artifact.get("order") < ORDER:
        raise ValueError("exact series is shorter than --order")
    exact_centers = {
        "D": QQ(-1)/2,
        "P": QQ(9)/4,
        "Q": -QQ(9)/4,
        "E": -QQ(27)/32,
    }
    exact_branch = {}
    for name in block_order:
        values = [QQ(value) for value in exact_artifact["series"][name]][:ORDER]
        values[0] -= exact_centers[name]
        exact_branch[name] = values
    assert exact_branch["P"] == [QQ(0), QQ(1)]+[QQ(0)]*(ORDER-2)
    for values in exact_branch.values():
        for value in values:
            if value.denominator() % 7 == 0:
                raise ValueError("exact branch series is not 7-integral")
    for name in block_order:
        reduced = [
            finite(value.numerator())/finite(value.denominator())
            for value in exact_branch[name]
        ]
        if reduced != branch_series[name]:
            disagreements = [
                index for index, (left, right) in enumerate(
                    zip(reduced, branch_series[name])
                ) if left != right
            ]
            raise ValueError(
                f"exact {name} jet does not reduce to the parameter branch "
                f"at coefficients {disagreements}"
            )

    # A fixed nonsingular minor is a left inverse for each base-7 Hensel digit.
    hensel_rows = tuple(jacobian.transpose().pivots())
    hensel_inverse = jacobian.matrix_from_rows(hensel_rows).inverse()
    assert len(hensel_rows) == len(seed)
    lifted = [ZZ(value) for value in seed]
    modulus = ZZ(7)
    reconstructed = None

    for digit in range(1, arguments.max_hensel_digits+1):
        next_modulus = modulus*7
        modular = Zmod(next_modulus)

        def reduce_exact(value):
            value = QQ(value)
            return (
                modular(value.numerator())
                /modular(value.denominator())
            )

        modular_branch = {
            name: [reduce_exact(value) for value in values]
            for name, values in exact_branch.items()
        }
        residual = equations_over(
            [modular(value) for value in lifted], modular_branch, modular
        )
        quotient = vector(
            finite,
            [
                finite((ZZ(value.lift())//modulus) % 7)
                for value in residual
            ],
        )
        delta = -hensel_inverse*vector(
            finite, [quotient[index] for index in hensel_rows]
        )
        if jacobian*delta != -quotient:
            augmented_rank = jacobian.augment((-quotient).column()).rank()
            inconsistent = tuple(
                index for index, value in enumerate(jacobian*delta+quotient)
                if value
            )
            print(
                "Q80SURFACEPARAMQQ|slope=8/87|"
                f"order={ORDER}|hensel_digit={digit}|"
                f"jacobian_rank={jacobian.rank()}|"
                f"augmented_rank={augmented_rank}|"
                f"inconsistent_equations={','.join(map(str, inconsistent))}|"
                "status=FAIL_FIRST_ORDER_HENSEL_OBSTRUCTION",
                flush=True,
            )
            raise ArithmeticError(
                "the rigid GF(7) parameter does not lift with these degrees "
                "and gauge"
            )
        lifted = [
            (value+modulus*ZZ(correction)) % next_modulus
            for value, correction in zip(lifted, delta)
        ]
        modulus = next_modulus

        if digit % 5 and digit != arguments.max_hensel_digits:
            continue
        try:
            candidate = [
                QQ(ZZ(value).rational_reconstruction(modulus))
                for value in lifted
            ]
        except ArithmeticError:
            continue
        exact_residual = equations_over(candidate, exact_branch, QQ)
        if not any(exact_residual):
            reconstructed = candidate
            break

    if reconstructed is None:
        raise ArithmeticError(
            f"no exact reconstruction after {arguments.max_hensel_digits} "
            "base-7 digits"
        )

    qq_polynomials = PolynomialRing(QQ, "t")
    qq_t = qq_polynomials.gen()
    parameter_functions_qq = {}
    for name in block_order:
        numerator, denominator = unpack(reconstructed, name)
        parameter_functions_qq[name] = (
            sum(value*qq_t**(index+1) for index, value in enumerate(numerator))
            /(
                1+sum(
                    value*qq_t**(index+1)
                    for index, value in enumerate(denominator)
                )
            )
        )

    # Re-expand the recovered functions at the CM point and compare all 28
    # exact branch coefficients, not only the independent Hensel rows.
    recovered_parameter_series = rational_series(
        [QQ(0)]+list(unpack(reconstructed, "P")[0]),
        [QQ(1)]+list(unpack(reconstructed, "P")[1]),
        ORDER,
    )
    for name in ("D", "Q", "E"):
        numerator, denominator = unpack(reconstructed, name)
        recovered = rational_series(
            [QQ(0)]+list(numerator),
            [QQ(1)]+list(denominator),
            ORDER,
        )
        expected = compose_series(
            exact_branch[name], recovered_parameter_series, ORDER
        )
        assert recovered == expected

    output = Path(arguments.qq_output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(
            {
                "schema": "q80-cm24-rational-surface-parameter-v1",
                "status": "exact_order_28_parameter_reconstruction",
                "scope": "candidate_global_parameter_pending_direct_surface_substitution",
                "slope": "8/87",
                "parameter_at_cm24": "t=0",
                "gauge": "P_num_t=-1,P_den_t=-2",
                "hensel_prime": 7,
                "hensel_digits": digit,
                "hensel_modulus_bits": modulus.nbits(),
                "source_series": str(series_path),
                "functions": {
                    name: {
                        "degrees": list(DEGREES[name]),
                        "value": str(value),
                    }
                    for name, value in parameter_functions_qq.items()
                },
            },
            indent=2,
            sort_keys=True,
        )+"\n"
    )
    largest_bits = max(
        max(abs(value.numerator()).nbits(), value.denominator().nbits())
        for value in reconstructed
    )
    print(
        "Q80SURFACEPARAMQQ|slope=8/87|"
        f"order={ORDER}|hensel_digits={digit}|modulus_bits={modulus.nbits()}|"
        f"largest_numden_bits={largest_bits}|output={output}|"
        "status=PASS_EXACT_PARAMETER_RECONSTRUCTION",
        flush=True,
    )
